<? // -*- html -*-
require_once("docutil.php");
page_head("Setting up a single-host server");
?>

<p>
BOINC provides a set of scripts for setting up
and controlling a BOINC server complex.
These scripts require all the server components to run on a single host.
This has been tested only on Linux and Solaris hosts;
it may work with small modifications on Windows also.
<p>
The scripts can be used to create multiple BOINC projects on the same host.
This can be handy for creating separate projects for testing and debugging.
In fact, the scripts are part of a general
<a href=test.html>testing framework</a>
that allows multiple developers to work independently on
a single host, with each developer able to create multiple projects.
<p>

Install all components listed in the <a href=software.php>Software
Prerequisites</a> page.  Your operating system must have shared memory
enabled, with a max shared segment size of at least 32 MB.

<h3>Creating the server</h3>
<p>
  Run the <code>make_project</code> script; example command line:
  <pre>
    cd tools/
    ./make_project --base $HOME/boinc --url_base http://boink/ yah 'YETI @ Home' upper_case 'UpperCase'
  </pre>

  See 'make_project --help' for more command-line options available (such as
  finer control of directory structure or clobbering an existing installation).
<p>
  The script does the following:
  <ul>
    <li> Create the project directory and its subdirectories.
    <li> Create the project's encryption keys if necessary.
      NOTE: before making the project publicly visible,
      you must move the code-signing private key
      to a highly secure (preferably non-networked) host,
      and delete it from the server host.
    <li> Create a MySQL database for this project,
      named (BOINC_USER_NAME)_(project-name).
    <li> Insert initial records in the project, platform, user, and app
         tables.
    <li> Copy cgi-programs to project/cgi-bin
    <li> Copy daemons to project/cgi
    <li> Copy html files to project/html_user and project/html_ops
    <li> Generate the configuration file (config.xml) used by
      the server programs.
  </ul>

<h3>Directory structure</h3>

Designate a 'BOINC projects' directory on the server host.
The scripts will create subdirectories as follows:
<pre>
boinc_projects/
    proj1/
        config.xml
        bin/
        cgi-bin/
        log/
        pid/
        download/
        html_ops/
        html_user/
        keys/
        upload/
    proj2/
    ...
</pre>
where proj1, proj2 etc. are the names of the projects you create.
Each project directory contains:
<ul>
  <li>config.xml: configuration file
  <li>bin: server daemons and programs, including the main <code>start</code>
    program as well as all the daemons
  <li>cgi-bin: cgi programs
  <li>log: log output
  <li>pid: lock files and pid files
  <li> download: storage for data server downloads.
  <li> html_ops: copies of PHP files for project management.
  <li> html_user: copies of PHP files for the public web site.
  <li> keys: encryption keys used by the project.
  <li> upload: storage for data server uploads.
</ul>

When you run the <code>make_project</code> script it will give you lines to
append to your Apache <code>httpd.conf</code>.  (Basically you need to alias
html_user, alias html_ops, and script-alias cgi-bin, all with appropriate
directory permissions.)  It will also give you a crontab line to install.

<p>
You should also set the default MIME type as follows:
<pre>
DefaultType application/octet-stream
</pre>


<p>
At this point you hopefully have a functioning BOINC server, but it has no
applications or work to distribute.  The remaining steps to make a public
project include:
<ul>
<li> Develop, debug and test your application.
<li> Develop back-end systems for generating work and processing results.
<li> Using the 'add' utility, add application versions to the BOINC database.
<li> Develop your web site.
</ul>

<?
page_tail();
?>
