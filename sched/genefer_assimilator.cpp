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

// A sample assimilator that:
// 1) if success, copy the output file(s) to a directory
// 2) if failure, append a message to an error log

#include <vector>
#include <string>
#include <cstdlib>

#include "boinc_db.h"
#include "error_numbers.h"
#include "filesys.h"
#include "sched_msgs.h"
#include "validate_util.h"
#include "sched_config.h"

using std::vector;
using std::string;

#define OUTPUT_DIR "genefer_results"

int write_error(char* p) {
    static FILE* f = 0;
    if (!f) {
        f = fopen(config.project_path(OUTPUT_DIR "/errors"), "a");
        if (!f) return ERR_FOPEN;
    }
    fprintf(f, "%s", p);
    fflush(f);
    return 0;
}

int get_residue_strdup(RESULT &result, const char *infile, char* &data, unsigned & EXP, unsigned & N);

int assimilate_handler(WORKUNIT& wu, vector<RESULT>& results, RESULT& canonical_result)
{
    int retval;
    char buf[1024];

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
                "[RESULT#%lu %s] unexpected number of input files (%d)\n",
                canonical_result.id, canonical_result.name, nFiles
            );
            return ERR_READ;
        }

        OUTPUT_FILE_INFO& fi = output_files[0];
        char *res = NULL;
        unsigned N, EXP;
        retval = get_residue_strdup(canonical_result, fi.path.c_str(), res, EXP, N);
        if (retval)
            return retval;
        unsigned EXP_SCALE;
        switch (EXP)
        {
        case  8192: EXP_SCALE = 13; break;
        case 16384: EXP_SCALE = 14; break;
        default:    EXP_SCALE = 99; break;
        }
        time_t now = time(NULL);
        strftime(buf, sizeof(buf), "%F", gmtime(&now));
        const char *output_path = config.project_path(OUTPUT_DIR "/gfn%u_%s.csv", EXP_SCALE, buf);
        FILE *f = fopen(output_path, "at");
        if (f == NULL)
        {
            log_messages.printf(MSG_CRITICAL,
                "[RESULT#%lu %s] cannot open output file '%s'\n",
                canonical_result.id, canonical_result.name, output_path
            );
            return ERR_FOPEN;
        }
        bool err = false;
        if (fprintf(f, "%u,%u,%s", EXP, N, res) < 0)
            err = true;
        for (unsigned i=0; i < results.size(); i++)
        {
            RESULT &r = results[i];
            if (r.validate_state == VALIDATE_STATE_VALID)
            {
                if (fprintf(f, ",%ld,%ld,%ld,%d", r.userid, r.hostid, r.app_version_id, r.received_time) < 0)
                    err = true;
            }
        }
        if (fprintf(f, "\n") < 0)
            err = true;

        if (fclose(f) < 0)
            err = true;
        if (err)
        {
            log_messages.printf(MSG_CRITICAL,
                "[RESULT#%lu %s] cannot write to output file '%s'\n",
                canonical_result.id, canonical_result.name, output_path
            );
            return ERR_FWRITE;
        }
    }
    else
    {
        sprintf(buf, "%s: 0x%x\n", wu.name, wu.error_mask);
        return write_error(buf);
    }
    return 0;
}
