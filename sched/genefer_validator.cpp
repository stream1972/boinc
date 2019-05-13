#include "config.h"
#include "util.h"
#include "sched_util.h"
#include "sched_msgs.h"
#include "validate_util.h"
#include "validate_util2.h"
#include "validator.h"
#include <inttypes.h>

using std::string;
using std::vector;

bool use_llr;

#define GFN16MEGA_B_OFFSET           1814570322693370

int validate_handler_init(int argc, char **argv)
{
    int i;

    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--llr") == 0)
            use_llr = true;
        else
        {
            fprintf(stderr, "Unknown option '%s'\n", argv[i]);
            return -1;
        }
    }
    return 0;
}

void validate_handler_usage(void)
{
    // describe the project specific arguments here
    fprintf(stderr,
        "    Custom options:\n"
        "    [--llr]                    parse LLR output files (default: Genefer)\n"
    );
}

static int get_residue_from_fp(RESULT &result, const char *infile, char* &data, unsigned &RET_EXP, uint64_t &RET_B, FILE *in)
{
    static bool once;
    static regex_t reg_composite, reg_prime;
    static regex_t reg_compo_llr, reg_prime_llr;
    if (!once)
    {
        // Genefer:
        // 13100008^8192+1 is composite. (RES=079b7b164cd29f71) (58305 digits) (err = 0.0000) (time = 0:00:09) 00:56:22
        // 13199802^8192+1 is a probable prime. (58332 digits) (err = 0.0044) (time = 0:01:32) 22:04:35
        // LLR:
        // 1814570322703210^65536+1 is not prime.  RES64: 34256CFC2ED1BF8E.  OLD64: 9C7046F48C753EA7  Time : 10446.820 sec.
        // 143946294^8192+1 is prime! (66832 decimal digits)  Time : 376.664 sec.
        static const char 
            pattern_composite[] = "^([0-9]+)\\^([0-9]+)\\+1 is composite\\. \\(RES=([0-9a-fA-F]{16})\\) \\([0-9]+ digits\\) \\(err = [0-9]",
            pattern_compo_llr[] = "^([0-9]+)\\^([0-9]+)\\+1 is not prime\\.  RES64: ([0-9a-fA-F]{16})\\.  OLD64: ([0-9a-fA-F]{16})  Time : [0-9]+",
            pattern_prime_llr[] =" ^([0-9]+)\\^([0-9]+)\\+1 is prime! \\([0-9]+ decimal digits\\)  Time : [0-9]+",
            pattern_prime[]     = "^([0-9]+)\\^([0-9]+)\\+1 is a probable prime\\. \\([0-9]+ digits\\) \\(err = [0-9]";

        once = true;
        if (regcomp(&reg_composite, pattern_composite, REG_EXTENDED | REG_ICASE) ||
            regcomp(&reg_compo_llr, pattern_compo_llr, REG_EXTENDED | REG_ICASE) ||
            regcomp(&reg_prime_llr, pattern_prime_llr, REG_EXTENDED | REG_ICASE) ||
            regcomp(&reg_prime,     pattern_prime,     REG_EXTENDED | REG_ICASE)
           )
        {
            log_messages.printf(MSG_CRITICAL,
                "[RESULT#%lu %s] regex compile failed\n",
                result.id, result.name
            );
            exit(1);
        }
    }

    int retval = ERR_OPENDIR;   // transient error to diagnose later
    char inbuf[1024];

read_again:
    if (fgets(inbuf, sizeof(inbuf), in) == NULL)
    {
        if (data)
            return 0;

        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] no signatures found in input file '%s'\n",
            result.id, result.name, infile
        );
        return retval;
    }

    char *p;
    while ((p = strpbrk(inbuf, "\r\n")) != NULL)  // CR-LF to spaces
        *p = ' ';
    for (int sz=strlen(inbuf); sz;)  // remove trailing spaces
    {
        if (inbuf[--sz] != ' ') break;
        inbuf[sz] = 0;
    }

    #define PMATCH_COUNT 5
    bool prime;
    regmatch_t pmatch[PMATCH_COUNT];

    if (regexec((use_llr ? &reg_compo_llr : &reg_composite), inbuf, PMATCH_COUNT, pmatch, 0) == 0)
    {
        prime = false;
    }
    else if (regexec((use_llr ? &reg_prime_llr : &reg_prime), inbuf, PMATCH_COUNT, pmatch, 0) == 0)
    {
        prime = true;
    }
    else
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] input line not parsed: '%s'\n",
            result.id, result.name, inbuf
        );
        goto read_again;
    }

    // 0 - whole string
    // 1 - B
    // 2 - exponent
    // 3 - residue if not prime

    // N-check. guess it from result name
    const char *pp;
    int file_exp;
    pp = strrchr(infile, '/');
    pp = pp ? pp + 1 : infile;   // strip path if present
    if (memcmp(pp, "gfn", 3) || !isdigit(pp[3]))
        goto fail_name;
    file_exp = atoi(pp+3);
    pp = strchr(pp, '_');
    if (pp && memcmp(pp, "_mega_", 6) == 0)  // gfn16_mega_
        pp += 5;
    if (!pp || !isdigit(pp[1]))
    {
    fail_name:
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] cannot guess B from input file name '%s'\n",
            result.id, result.name, infile
        );
        return retval;
    }
    int64_t expected_value = atoll(pp+1);
    int64_t n = atoll(inbuf + pmatch[1].rm_so);
    if (file_exp == 16)
        expected_value += GFN16MEGA_B_OFFSET;
    if (n != expected_value)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] bad B %" PRIi64 ", expected %" PRIi64 "\n",
            result.id, result.name, n, expected_value
        );
        return retval;
    }
    RET_B = n;

    // exponent check
    switch (file_exp)
    {
    case 13: expected_value = 8192;   break;
    case 14: expected_value = 16384;  break;
    case 16: expected_value = 65536;  break;
    default:
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] cannot guess exponent from N %d\n",
            result.id, result.name, file_exp
        );
        return retval;
    }
    n = atoi(inbuf + pmatch[2].rm_so);
    if (n != expected_value)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] bad exponent %" PRIi64 ", expected %" PRIi64 "\n",
            result.id, result.name, n, expected_value
        );
        return retval;
    }
    RET_EXP = n;

    const char *curdata;
    if (prime)
        curdata = "PRIME";
    else
    {
        inbuf[pmatch[3].rm_eo] = 0;
        curdata = inbuf + pmatch[3].rm_so;
    }
    if (data == NULL)
    {
        data = strdup(curdata);
        for (p = data; *p; p++)
            *p = toupper(*p);
    }
    else if (strcmp(data, curdata))
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] multiple mismatched residues '%s' and '%s'\n",
            result.id, result.name,
            data, curdata
        );
        return retval;
    }
    goto read_again;
}

int get_residue_strdup(RESULT &result, const char *infile, char* &data, unsigned &RET_EXP, uint64_t &RET_B)
{
    FILE *in;

    int retval;

    data = NULL;
    retval = try_fopen(infile, in, "rt");
    if (retval)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] cannot open file '%s'\n",
            result.id, result.name, infile
        );
    }
    else
    {
        retval = get_residue_from_fp(result, infile, data, RET_EXP, RET_B, in);
        fclose(in);
    }

    return retval;
}

unsigned get_llr_version(const char *text)
{
    static const char pattern[] = "\nLLR Program - Version ([0-9]+)\\.([0-9]+)\\.([0-9]+),";
    static bool once;
    static regex_t reg;

    if (!once)
    {
        once = true;
        if (regcomp(&reg, pattern, REG_EXTENDED | REG_ICASE))
        {
            log_messages.printf(MSG_CRITICAL,
                "%s: regex compile failed\n", __FUNCTION__
            );
            exit(1);
        }
    }

    regmatch_t pmatch[PMATCH_COUNT];
    unsigned version = 0;
    if (regexec(&reg, text, PMATCH_COUNT, pmatch, 0) == 0)
    {
        int i;
        for (i = 1; i <= 3; i++)
            version = version * 100 + atoi(text + pmatch[i].rm_so);
    }
    return version;
}

int init_result(RESULT& result, void*& data)
{
    int retval;
    vector<OUTPUT_FILE_INFO> files;

    retval = get_output_file_infos(result, files);
    if (retval)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] check_set: can't get output filenames\n",
            result.id, result.name
        );
        return retval;
    }

    int nFiles = files.size();
    if (nFiles != 1)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] unexpected number of input files (%d)\n",
            result.id, result.name, nFiles
        );
        return ERR_OPENDIR;   // transient error to diagnose later
    }

    for (int i = 0; i < nFiles; i++)
    {
        OUTPUT_FILE_INFO& fi = files[i];

        if (fi.no_validate)
        {
            log_messages.printf(MSG_CRITICAL,
                "[RESULT#%lu %s] unexpected no_validate flag\n",
                result.id, result.name
            );
            return ERR_OPENDIR;   // transient error to diagnose later
        }

        char *res = NULL;
        unsigned dummy_exp;
        uint64_t dummy_b;
        retval = get_residue_strdup(result, fi.path.c_str(), res, dummy_exp, dummy_b);
//      printf("Got residue: '%s'\n", res);
        char *full_res = NULL;
        if (res)
        {
            full_res = (char*) malloc(strlen(res) + 64);
            if (full_res)
                sprintf(full_res, "%u %" PRIu64 " %s", dummy_exp, dummy_b, res);
            free(res);
        }
//      printf("Got data: '%s'\n", full_res);
        data = full_res;
        return retval;
    }

    return ERR_OPENDIR;
}

int compare_results(
    RESULT       & r1, void* data1,
    RESULT const & r2, void* data2,
    bool& match
    )
{
    if (data1 == NULL || data2 == NULL)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] NULL data in compare_results\n",
            r1.id, r1.name
        );
        return ERR_OPENDIR;   // transient error to diagnose later
    }

    match = strcasecmp((char*)data1, (char*)data2) == 0;
    if (!match)
    {
        log_messages.printf(MSG_CRITICAL, "Validation error:\n");
        log_messages.printf(MSG_CRITICAL, "  [RESULT#%lu %s U%-4ld H%-4ld] %s\n", r1.id, r1.name, r1.userid, r1.hostid, (char*)data1);
        log_messages.printf(MSG_CRITICAL, "  [RESULT#%lu %s U%-4ld H%-4ld] %s\n", r2.id, r2.name, r2.userid, r2.hostid, (char*)data2);
    }

    return 0;
}

int cleanup_result(RESULT const& /*result*/, void* data)
{
    if (data)
        free(data);
    return 0;
}
