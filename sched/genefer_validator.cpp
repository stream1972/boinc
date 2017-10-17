#include "config.h"
#include "util.h"
#include "sched_util.h"
#include "sched_msgs.h"
#include "validate_util.h"
#include "validate_util2.h"
#include "validator.h"

using std::string;
using std::vector;

static void parse_cmdline(void)
{
    for (int i=1; i<g_argc; i++)
    {
        // check g_argv[i] here
    }
}

static int get_residue_from_fp(RESULT &result, const char *infile, char* &data, unsigned &RET_EXP, unsigned &RET_N, FILE *in)
{
    static bool once;
    static regex_t reg_composite, reg_prime;
    if (!once)
    {
        // 13100008^8192+1 is composite. (RES=079b7b164cd29f71) (58305 digits) (err = 0.0000) (time = 0:00:09) 00:56:22
        // 13199802^8192+1 is a probable prime. (58332 digits) (err = 0.0044) (time = 0:01:32) 22:04:35
        static const char 
            pattern_composite[] = "^([0-9]+)\\^([0-9]+)\\+1 is composite\\. \\(RES=([0-9a-fA-F]{16})\\) \\([0-9]+ digits\\) \\(err = [0-9]",
            pattern_prime[]     = "^([0-9]+)\\^([0-9]+)\\+1 is a probable prime\\. \\([0-9]+ digits\\) \\(err = [0-9]";

        once = true;
        if (regcomp(&reg_composite, pattern_composite, REG_EXTENDED | REG_ICASE) ||
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
    int size;
    char inbuf[1024];
    if ((size = fread(inbuf, 1, sizeof(inbuf), in)) <= 0)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] error reading input file '%s'\n",
            result.id, result.name, infile
        );
        return retval;
    }
    if (size >= (int)sizeof(inbuf))
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] input file '%s' is too big\n",
            result.id, result.name, infile
        );
        return retval;
    }

    {
    char *p;
    inbuf[size] = 0;
    while ((p = strpbrk(inbuf, "\r\n")) != NULL)  // CR-LF to spaces
        *p = ' ';
    }

    #define PMATCH_COUNT 5
    bool prime;
    regmatch_t pmatch[PMATCH_COUNT];

    if (regexec(&reg_composite, inbuf, PMATCH_COUNT, pmatch, 0) == 0)
    {
        prime = false;
    }
    else if (regexec(&reg_prime, inbuf, PMATCH_COUNT, pmatch, 0) == 0)
    {
        prime = true;
    }
    else
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] input data not parsed: '%s'\n",
            result.id, result.name, inbuf
        );
        return retval;
    }

    // 0 - whole string
    // 1 - N
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
    if (!pp || !isdigit(pp[1]))
    {
    fail_name:
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] cannot guess N from input file name '%s'\n",
            result.id, result.name, infile
        );
        return retval;
    }
    int expected_value = atoi(pp+1);
    int n = atoi(inbuf + pmatch[1].rm_so);
    if (n != expected_value)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] bad N %d, expected %d\n",
            result.id, result.name, n, expected_value
        );
        return retval;
    }
    RET_N = n;

    // exponent check
    switch (file_exp)
    {
    case 13: expected_value = 8192;   break;
    case 14: expected_value = 16384;  break;
    default:
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] cannot guess exponent from N %d\n",
            result.id, result.name, n
        );
        return retval;
    }
    n = atoi(inbuf + pmatch[2].rm_so);
    if (n != expected_value)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] bad exponent %d, expected %d\n",
            result.id, result.name, n, expected_value
        );
        return retval;
    }
    RET_EXP = n;

    if (prime)
    {
        data = strdup("PRIME");
    }
    else
    {
        inbuf[pmatch[3].rm_eo] = 0;
        data = strdup(inbuf + pmatch[3].rm_so);
        for (char *p = data; *p; p++)
            *p = toupper(*p);
    }
    return 0;
}

int get_residue_strdup(RESULT &result, const char *infile, char* &data, unsigned &RET_EXP, unsigned &RET_N)
{
    FILE *in;

    int retval;

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
        retval = get_residue_from_fp(result, infile, data, RET_EXP, RET_N, in);
        fclose(in);
    }

    return retval;
}

int init_result(RESULT& result, void*& data)
{
    int retval;
    vector<OUTPUT_FILE_INFO> files;

    static bool once;
    if (!once)
    {
        parse_cmdline();
        once = true;
    }

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
        unsigned dummy_exp, dummy_n;
        retval = get_residue_strdup(result, fi.path.c_str(), res, dummy_exp, dummy_n);
//      printf("Got residue: '%s'\n", res);
        char *full_res = NULL;
        if (res)
            full_res = (char*) malloc(strlen(res) + 64);
        if (full_res)
            sprintf(full_res, "%u %u %s", dummy_exp, dummy_n, res);
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
    (void) r2;
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
        log_messages.printf(MSG_CRITICAL, "  [RESULT#%lu %s] %s\n", r1.id, r1.name, (char*)data1);
        log_messages.printf(MSG_CRITICAL, "  [RESULT#%lu %s] %s\n", r2.id, r2.name, (char*)data2);
    }

    return 0;
}

int cleanup_result(RESULT const& /*result*/, void* data)
{
    if (data)
        free(data);
    return 0;
}
