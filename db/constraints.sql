
alter table platform
    add unique(name);

alter table app
    add unique(name);

alter table app_version
    add unique apvp (appid, platformid, version_num, plan_class);

alter table user
    add unique(email_addr),
    add unique(authenticator),
    add index ind_tid (teamid),
    add index user_name(name),
    add index user_tot (total_credit desc),
        -- db_dump.C
    add index user_avg (expavg_credit desc);
        -- db_dump.C

alter table team
    add unique(name),
    add fulltext index team_name_desc(name, description),
    add fulltext index team_name(name),
    add index team_avg (expavg_credit desc),
        -- db_dump.C
    add index team_tot (total_credit desc),
        -- db_dump.C
    add index team_userid (userid);

alter table workunit
    add unique(name),
        -- not currently used but good invariant
    add index wu_val (appid, need_validate),
        -- validator
    add index wu_timeout (transition_time),
        -- transitioner
    add index wu_filedel (file_delete_state),
        -- file_deleter, db_purge
    add index wu_assim (appid, assimilate_state);
        -- assimilator

alter table result
    add unique(name),
        -- the scheduler looks up results by name

    add index res_wuid (workunitid),
        -- transitioner
        -- NOTE res_wu_user may suffice, could try dropping this one

    add index ind_res_st (server_state, priority),
        -- feeder

    add index res_app_state(appid, server_state),
        -- to get count of unsent results for given app (e.g. in work generator)

    add index res_filedel (file_delete_state),
        -- file_deleter

    add index res_userid_id(userid, id desc),
        -- html_user/results.php

    add index res_userid_val(userid, validate_state),
        -- to show pending credit

    add index res_hostid_id (hostid, id desc),
        -- html_user/results.php

    add index res_wu_user (workunitid, userid),
        -- scheduler (avoid sending mult results of same WU to one user)

    -- my extra fields

    add index res_senttime (sent_time),
        -- HTML: ordering of workunits (???)

    add index res_user_gui     (userid, gui_state),
    add index res_user_app_gui (userid, appid, gui_state),
    add index res_host_gui     (hostid, gui_state),
    add index res_host_app_gui (hostid, appid, gui_state);
        -- HTML: quick calc for number of results for user/app/host combos

alter table msg_from_host
    add index message_handled (handled);
        -- for message handler

alter table msg_to_host
    add index msg_to_host(hostid, handled);
        -- for scheduler

alter table host
    add index host_user (userid),
        -- html_user/host_user.php
    add index host_avg (expavg_credit desc),
        -- db_dump.C
    add index host_tot (total_credit desc);
        -- db_dump.C

alter table profile
    add fulltext index profile_reponse(response1, response2),
    add index pro_uotd (uotd_time desc),
    add unique profile_userid(userid);

alter table subscriptions
    add unique sub_unique(userid, threadid);

alter table category
    add unique cat1(name, is_helpdesk);

alter table forum
    add unique pct (parent_type, category, title);

alter table thread
    add fulltext index thread_title(title);
        
alter table post
    add index post_user (user),
    add index post_thread (thread),
    add fulltext index post_content(content);

alter table credited_job 
    add index credited_job_user (userid),
    add index credited_job_wu (workunitid),
    add unique credited_job_user_wu (userid, workunitid);

alter table team_delta
    add index team_delta_teamid (teamid, timestamp);

alter table team_admin
    add unique (teamid, userid);

alter table friend
    add unique friend_u (user_src, user_dest);

alter table notify
    add unique notify_un (userid, type, opaque);

alter table host_app_version
    add unique hap(host_id, app_version_id);

alter table assignment
    add index asgn_target(target_type, target_id);

alter table job_file
    add unique jf_name(name);

alter table badge_user
    add unique (user_id, badge_id);

alter table badge_team
    add unique (team_id, badge_id);

alter table credit_user
    add index cu_total(appid, total),
    add index cu_avg(appid, expavg);

alter table credit_team
    add index ct_total(appid, total),
    add index ct_avg(appid, expavg);




    -- gui_state field must be calculated by database itself
delimiter //

create trigger upd_gui_state before update on result
for each row
begin
    if NEW.server_state <> OLD.server_state OR NEW.outcome <> OLD.outcome OR NEW.validate_state <> OLD.validate_state then
        set NEW.gui_state =
            case NEW.server_state
                when 4 then 1
                when 5 then
                    case
                        when NEW.outcome=1 then
                            case NEW.validate_state
                                when 0 then 2
                                when 4 then 3
                                when 1 then 4
                                when 2 then 5
                                when 3 then 5
                                when 5 then 5
                                else NEW.gui_state
                            end
                        when NEW.outcome=6            then 5
                        when NEW.outcome in (3, 4, 7) then 6
                        else NEW.gui_state
                    end
                else NEW.gui_state
            end
        ;
    end if;
end; //

delimiter ;
