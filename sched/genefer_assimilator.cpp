// This file is part of BOINC.
// http://boinc.berkeley.edu
// Copyright (C) 2008 University of California
//
// BOINC is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// BOINC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with BOINC.  If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <string>
#include <cstdlib>

#include "boinc_db.h"
#include "error_numbers.h"
#include "filesys.h"
#include "sched_msgs.h"
#include "validate_util.h"
#include "sched_config.h"
#include <inttypes.h>

using std::vector;
using std::string;

#define OUTPUT_DIR "genefer_results"

int get_residue_strdup(RESULT &result, const char *infile, char* &data, unsigned & EXP, uint64_t & B);
unsigned get_llr_version(const char *text);
extern bool use_llr;

int assimilate_handler_init(int argc, char** argv)
{
    // handle project specific arguments here
    return validate_handler_init(argc, argv);
}

void assimilate_handler_usage(void)
{
    // describe the project specific arguments here
    //fprintf(stderr,
    //    "    Custom options:\n"
    //    "    [--project_option X]  a project specific option\n"
    //);
    validate_handler_usage();
}

static int copy_result_files(RESULT &result, FILE *fout)
{
    char buf[1024];
    vector<OUTPUT_FILE_INFO> files;
    int retval;

    retval = get_output_file_infos(result, files);
    if (retval)
    {
        log_messages.printf(MSG_CRITICAL,
            "[RESULT#%lu %s] can't get output filenames\n",
            result.id, result.name
        );
        return retval;
    }

    int nFiles = files.size();
    for (int i = 0; i < nFiles; i++)
    {
        OUTPUT_FILE_INFO& fi = files[i];
        FILE *f = fopen(fi.path.c_str(), "rt");
        if (f == NULL)
        {
            if (result.validate_state != VALIDATE_STATE_VALID)
                continue;
            log_messages.printf(MSG_CRITICAL,
                "[RESULT#%lu %s] can't open '%s'%s\n",
                result.id, result.name,
                fi.path.c_str(),
                (result.validate_state == VALIDATE_STATE_VALID ? "" : ", ignored for invalid result")
                );
            return ERR_FOPEN;
        }
        if (fprintf(fout, "[%s]\n", fi.name.c_str()) < 0)
            retval = ERR_FWRITE;
        while (!retval && fgets(buf, sizeof(buf)-1, f))
        {
            int sz = strlen(buf);
            // force \n on corrupted lines
            if (sz && buf[sz-1] != '\n')
            {
                buf[sz] = '\n';
                buf[sz+1] = 0;
            }
            if (fputs(buf, fout) == EOF)
                retval = ERR_FWRITE;
        }
        fclose(f);
    }

    return retval;
}

int assimilate_handler(WORKUNIT& wu, vector<RESULT>& results, RESULT& canonical_result)
{
    int retval;
    char tsbuf[128];

    static bool once;
    if (!once)
    {
        once = true;
        retval = boinc_mkdir(config.project_path(OUTPUT_DIR));
        if (retval) return retval;
    }

    if (wu.canonical_resultid)
    {
        vector<OUTPUT_FILE_INFO> output_files;
        get_output_file_infos(canonical_result, output_files);
        int nFiles = output_files.size();
        if (nFiles != 1)
        {
            log_messages.printf(MSG_CRITICAL,
                "[RESULT#%lu %s] unexpected number of files (%d)\n",
                canonical_result.id, canonical_result.name, nFiles
            );
            return ERR_READ;
        }

        OUTPUT_FILE_INFO& fi = output_files[0];
        char *res = NULL;
        unsigned EXP;
        uint64_t B;
        retval = get_residue_strdup(canonical_result, fi.path.c_str(), res, EXP, B);
        if (retval)
            return retval;
        unsigned EXP_SCALE;
        switch (EXP)
        {
        case  8192: EXP_SCALE = 13; break;
        case 16384: EXP_SCALE = 14; break;
        case 65536: EXP_SCALE = 16; break;
        default:    EXP_SCALE = 99; break;
        }
        time_t now = time(NULL);
        strftime(tsbuf, sizeof(tsbuf), (wu.transitioner_flags & TRANSITION_NO_NEW_RESULTS ? "%Y-%m-cancelled" : "%F"), gmtime(&now));
        const char *output_path = config.project_path(OUTPUT_DIR "/gfn%u_%s.csv%s", EXP_SCALE, tsbuf, (use_llr ? "_new" : ""));
        FILE *f = fopen(output_path, "at");
        if (f == NULL)
        {
            log_messages.printf(MSG_CRITICAL,
                "[RESULT#%lu %s] cannot open output file '%s'\n",
                canonical_result.id, canonical_result.name, output_path
            );
            free(res);
            return ERR_FOPEN;
        }

        bool err = false;
        if (fprintf(f, "%u,%" PRIu64 ",%s", EXP, B, res) < 0)
            err = true;
        for (unsigned i=0; i < results.size(); i++)
        {
            RESULT &r = results[i];
            if (r.validate_state == VALIDATE_STATE_VALID)
            {
                if (use_llr)
                {
                    if (fprintf(f, ",%ld,%ld,%ld,%ld,%d,%.15g", r.userid, r.hostid, r.teamid, r.app_version_id, r.received_time, r.granted_credit) < 0)
                        err = true;
                    if (fprintf(f, ",%u", get_llr_version(r.stderr_out)) < 0)
                        err = true;
                }
                else
                {
                    if (fprintf(f, ",%ld,%ld,%ld,%d", r.userid, r.hostid, r.app_version_id, r.received_time) < 0)
                        err = true;
                }
            }
        }
        if (fprintf(f, "\n") < 0)
            err = true;

        if (fclose(f) < 0)
            err = true;

        free(res);
        if (err)
        {
            log_messages.printf(MSG_CRITICAL,
                "[RESULT#%lu %s] cannot write to output file '%s'\n",
                canonical_result.id, canonical_result.name, output_path
            );
            return ERR_FWRITE;
        }

        // Keep extra archive of full LLR results
        if (use_llr)
        {
            FILE *fout;

            output_path = config.project_path(OUTPUT_DIR "/residues-gfn%u_%s.txt", EXP_SCALE, tsbuf);
            fout = fopen(output_path, "at");
            if (fout == NULL)
            {
                log_messages.printf(MSG_CRITICAL,
                    "[RESULT#%lu %s] cannot open output file '%s'\n",
                    canonical_result.id, canonical_result.name, output_path
                );
                return ERR_FOPEN;
            }

            retval = 0;
            for (unsigned i = 0; i < results.size() && !retval; i++)
                retval = copy_result_files(results[i], fout);

            if (fclose(fout) < 0)
                retval = ERR_FWRITE;

            if (retval)
            {
                log_messages.printf(MSG_CRITICAL,
                    "[RESULT#%lu %s] cannot write to output file '%s'\n",
                    canonical_result.id, canonical_result.name, output_path
                );
                return retval;
            }
        }
    }
    else
    {
        const char *output_path = config.project_path(OUTPUT_DIR "/errors");
        FILE *f = fopen(output_path, "at");
        if (f == NULL)
        {
            log_messages.printf(MSG_CRITICAL, "cannot open output file '%s'\n", output_path);
            return ERR_FOPEN;
        }
        bool err = false;
        if (fprintf(f, "%s: 0x%x\n", wu.name, wu.error_mask) < 0)
            err = true;
        if (fclose(f) < 0)
            err = true;
        if (err)
        {
            log_messages.printf(MSG_CRITICAL, "cannot write to output file '%s'\n", output_path);
            return ERR_FWRITE;
        }
    }
    return 0;
}
