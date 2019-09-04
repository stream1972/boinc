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

// job dispatch using a score-based approach:
// - scan the job array, assigning a score to each job and building a list
//   (the score reflect a variety of factors).
// - sort the list
// - send jobs in order of decreasing score until request is satisfied
// - do the above separately for each resource type

#include <algorithm>

#include "boinc_db.h"
#include "error_numbers.h"
#include "util.h"

#include "sched_check.h"
#include "sched_config.h"
#include "sched_hr.h"
#include "sched_keyword.h"
#include "sched_main.h"
#include "sched_msgs.h"
#include "sched_send.h"
#include "sched_shmem.h"
#include "sched_types.h"
#include "sched_version.h"

#include "sched_score.h"

// given the host's estimated speed, determine its size class
//
static int get_size_class(APP& app, double es) {
    for (int i=0; i<app.n_size_classes-1; i++) {
        if (es < app.size_class_quantiles[i]) return i;
    }
    return app.n_size_classes - 1;
}

// Assign a score to this job,
// representing the value of sending the job to this host.
// Also do some initial screening,
// and return false if can't send the job to host
//
bool JOB::get_score(int array_index) {
    WU_RESULT& wu_result = ssp->wu_results[array_index];
    score = 0;

    if (!app->beta && wu_result.need_reliable) {
        if (!bavp->reliable) {
            return false;
        }
    }

    if (app->beta) {
        if (g_wreq->project_prefs.allow_beta_work) {
            score += 1;
        } else {
            if (config.debug_send_job) {
                log_messages.printf(MSG_NORMAL,
                    "[send_job] can't send job %lu for beta app to non-beta user\n",
                    wu_result.workunit.id
                );
            }
            return false;
        }
    }

    if (app_not_selected(app->id)) {
        if (g_wreq->project_prefs.allow_non_preferred_apps) {
            score -= 1;
        } else {
            if (config.debug_send_job) {
                log_messages.printf(MSG_NORMAL,
                    "[send_job] app not selected for job %lu\n",
                    wu_result.workunit.id
                );
            }
            return false;
        }
    }

    // Skip versions (CPU/GPU) disabled on site in project-specific settings.
    int pt = bavp->host_usage.proc_type;
    if (proc_type_not_selected(app->id, pt)) {
        if (g_wreq->project_prefs.allow_non_preferred_apps) {
            score -= 1;
        } else {
            if (config.debug_send_job) {
                log_messages.printf(MSG_NORMAL,
                    "[send_job] proctype %d (%s) not selected for job %lu\n",
                    pt, proc_type_name(pt), wu_result.workunit.id
                );
            }
            return false;
        }
    }

    if (wu_result.infeasible_count) {
        score += 1;
    }

    if (app->locality_scheduling == LOCALITY_SCHED_LITE
        && g_request->file_infos.size()
    ) {
        int n = nfiles_on_host(wu_result.workunit);
        if (config.debug_locality_lite) {
            log_messages.printf(MSG_NORMAL,
                "[loc_lite] job %s has %d files on this host\n",
                wu_result.workunit.name, n
            );
        }
        if (n > 0) {
            score += 10;
        }
    }

    if (app->n_size_classes > 1) {
        double effective_speed = bavp->host_usage.projected_flops * available_frac(*bavp);
        int target_size = get_size_class(*app, effective_speed);
        if (config.debug_send_job) {
            log_messages.printf(MSG_NORMAL,
                "[send_job] size: host %d job %d speed %f\n",
                target_size, wu_result.workunit.size_class, effective_speed
            );
        }
        if (target_size == wu_result.workunit.size_class) {
            score += 5;
        } else if (target_size < wu_result.workunit.size_class) {
            score -= 2;
        } else {
            score -= 1;
        }
    }

#if 1
    // my way: prefer old workunits
    double d = (double)(time(NULL) - wu_result.workunit.create_time) / (365*24*60*60);
    if (d > 0)
        score += d;
#else
    // Boinc default way. Unacceptable for me, becase priority contains 'n' of candidate
    score += wu_result.res_priority;
#endif

    if (config.keyword_sched) {
        double x = keyword_score(array_index);
        if (x < 0) {
            return false;
        }
        score += x;
    }

    if (config.debug_send_job) {
        log_messages.printf(MSG_NORMAL,
            "[send_job]: score %f for result %lu\n", score, wu_result.resultid
        );
    }

    return true;
}

bool job_compare(JOB j1, JOB j2) {
    return (j1.score > j2.score);
}

static double req_sec_save[NPROC_TYPES];
static double req_inst_save[NPROC_TYPES];

static void clear_others(int rt) {
    for (int i=0; i<NPROC_TYPES; i++) {
        if (i == rt) continue;
        req_sec_save[i] = g_wreq->req_secs[i];
        g_wreq->req_secs[i] = 0;
        req_inst_save[i] = g_wreq->req_instances[i];
        g_wreq->req_instances[i] = 0;
    }
}

static void restore_others(int rt) {
    for (int i=0; i<NPROC_TYPES; i++) {
        if (i == rt) continue;
        g_wreq->req_secs[i] += req_sec_save[i];
        g_wreq->req_instances[i] += req_inst_save[i];
    }
}

// send work for a particular processor type
//
void send_work_score_type(int rt) {
    vector<JOB> jobs;

    if (config.debug_send_scan) {
        log_messages.printf(MSG_NORMAL,
            "[send_scan] scanning for %s jobs\n", proc_type_name(rt)
        );
    }

    clear_others(rt);

    int nscan = ssp->max_wu_results;
    int rnd_off = rand() % ssp->max_wu_results;
    if (config.debug_send_scan) {
        log_messages.printf(MSG_NORMAL,
            "[send_scan] scanning %d slots starting at %d\n", nscan, rnd_off
        );
    }
    for (int j=0; j<nscan; j++) {
        int i = (j+rnd_off) % ssp->max_wu_results;
        WU_RESULT& wu_result = ssp->wu_results[i];
        if (wu_result.state != WR_STATE_PRESENT  && wu_result.state != g_pid) {
            continue;
        }
        WORKUNIT wu = wu_result.workunit;
        JOB job;
        job.app = ssp->lookup_app(wu.appid);
        if (job.app->non_cpu_intensive) {
            if (config.debug_send_job) {
                log_messages.printf(MSG_NORMAL,
                    "[send_job] [RESULT#%lu] app is non compute intensive\n",
                    wu_result.resultid
                );
            }
            continue;
        }
        job.bavp = get_app_version(wu, true, false);
        if (!job.bavp) {
            if (config.debug_send_job) {
                log_messages.printf(MSG_NORMAL,
                    "[send_job] [RESULT#%lu] no app version available\n",
                    wu_result.resultid
                );
            }
            continue;
        }

        job.index = i;
        job.result_id = wu_result.resultid;
        if (!job.get_score(i)) {
            if (config.debug_send_job) {
                log_messages.printf(MSG_NORMAL,
                    "[send_job] [RESULT#%lu] get_score() returned false\n",
                    wu_result.resultid
                );
            }
            continue;
        }
        if (config.debug_send_job) {
            log_messages.printf(MSG_NORMAL,
                "[send_job] [RESULT#%lu] score: %f\n",
                wu_result.resultid, job.score
            );
        }
        jobs.push_back(job);
    }

    std::sort(jobs.begin(), jobs.end(), job_compare);

    bool sema_locked = false;
    for (unsigned int i=0; i<jobs.size(); i++) {

        // check limit on total jobs
        //
        if (!work_needed(false)) {
            break;
        }

        // check limit on jobs for this processor type
        //
        if (!g_wreq->need_proc_type(rt)) {
            break;
        }
        JOB& job = jobs[i];

        // check limits on jobs for this (app, processor type)
        //
        if (config.max_jobs_in_progress.exceeded(job.app, job.bavp->host_usage.proc_type)) {
            if (config.debug_quota) {
                log_messages.printf(MSG_NORMAL,
                    "[quota] limit for app/proctype exceeded\n"
                );
            }
            continue;
        }

        if (!sema_locked) {
            lock_sema();
            sema_locked = true;
        }

        // make sure the job is still in the cache
        // array is locked at this point.
        //
        WU_RESULT& wu_result = ssp->wu_results[job.index];
        if (wu_result.state != WR_STATE_PRESENT  && wu_result.state != g_pid) {
            continue;
        }
        if (wu_result.resultid != job.result_id) {
            continue;
        }
        WORKUNIT wu = wu_result.workunit;
        int retval = wu_is_infeasible_fast(
            wu,
            wu_result.res_server_state, wu_result.res_priority,
            wu_result.res_report_deadline,
            *job.app,
            *job.bavp
        );

        if (retval) {
            continue;
        }
        wu_result.state = g_pid;

        // It passed fast checks.
        // Release sema and do slow checks
        //
        unlock_sema();
        sema_locked = false;

        switch (slow_check(wu_result, job.app, job.bavp)) {
        case CHECK_NO_HOST:
            wu_result.state = WR_STATE_PRESENT;
            break;
        case CHECK_NO_ANY:
            wu_result.state = WR_STATE_EMPTY;
            if (config.keyword_sched) {
                keyword_sched_remove_job(job.index);
            }
            break;
        default:
            // slow_check() refreshes fields of wu_result.workunit;
            // update our copy too
            //
            wu.hr_class = wu_result.workunit.hr_class;
            wu.app_version_id = wu_result.workunit.app_version_id;

            // mark slot as empty AFTER we've copied out of it
            // (since otherwise feeder might overwrite it)
            //
            if (config.keyword_sched) {
                keyword_sched_remove_job(job.index);
            }

            // reread result from DB, make sure it's still unsent
            // TODO: from here to end of add_result_to_reply()
            // (which updates the DB record) should be a transaction
            //
            SCHED_DB_RESULT result;
            result.id = wu_result.resultid;
            if (result_still_sendable(result, wu)) {
                add_result_to_reply(result, wu, job.bavp, false);

                // add_result_to_reply() fails only in pathological cases -
                // e.g. we couldn't update the DB record or modify XML fields.
                // If this happens, don't replace the record in the array
                // (we can't anyway, since we marked the entry as "empty").
                // The feeder will eventually pick it up again,
                // and hopefully the problem won't happen twice.
            }

            // Free record only after DB update to decrease chance of race condition in feeder
            wu_result.state = WR_STATE_EMPTY;
            break;
        }
    }
    if (sema_locked) {
        unlock_sema();
    }

    restore_others(rt);
    g_wreq->best_app_versions.clear();
}

void send_work_score() {
    for (int i=NPROC_TYPES-1; i>= 0; i--) {
        if (g_wreq->need_proc_type(i)) {
            send_work_score_type(i);
        }
    }
}
