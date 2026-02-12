/*
 * GLINT Installer Launcher (macOS)
 *
 * Compiled Mach-O wrapper that replaces the bash launcher script.
 * When macOS launches an .app bundle whose CFBundleExecutable is a shell
 * script, Terminal.app opens as an intermediary and the tkinter GUI may
 * never reach the foreground. A native binary avoids this entirely.
 *
 * Build:
 *   clang -O2 -o launcher launcher.c -framework ApplicationServices
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <libgen.h>        /* dirname() */
#include <sys/stat.h>

#include <mach-o/dyld.h>   /* _NSGetExecutablePath */
#include <ApplicationServices/ApplicationServices.h>

/* ── Logging ─────────────────────────────────────────────────────── */

static FILE *g_log = NULL;

static void log_open(void) {
    const char *home = getenv("HOME");
    if (!home) home = "/tmp";

    char log_dir[1024];
    snprintf(log_dir, sizeof(log_dir), "%s/Library/Logs", home);
    mkdir(log_dir, 0755);                       /* ignore error if exists */

    char log_path[1024];
    snprintf(log_path, sizeof(log_path), "%s/GLINT_Installer.log", log_dir);

    g_log = fopen(log_path, "a");
}

static void log_msg(const char *fmt, ...) {
    if (!g_log) return;

    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    char ts[64];
    strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S", t);
    fprintf(g_log, "[%s] ", ts);

    va_list ap;
    va_start(ap, fmt);
    vfprintf(g_log, fmt, ap);
    va_end(ap);

    fprintf(g_log, "\n");
    fflush(g_log);
}

/* ── Utilities ───────────────────────────────────────────────────── */

/* Return 1 if `python_path` exists, is executable, can import tkinter,
 * and has Tk version >= 8.6.  Tk 8.5 (shipped with macOS system Python)
 * has rendering issues with ttk on modern macOS (Monterey+). */
static int python_has_modern_tk(const char *python_path) {
    if (access(python_path, X_OK) != 0)
        return 0;

    char cmd[1024];
    snprintf(cmd, sizeof(cmd),
             "%s -c 'import tkinter; exit(0 if float(tkinter.TkVersion) >= 8.6 else 1)'"
             " >/dev/null 2>&1", python_path);

    int rc = system(cmd);
    if (rc != 0)
        log_msg("    (Tk < 8.6 or tkinter missing)");
    return (rc == 0);
}

/* Weaker check: can import tkinter at all (even Tk 8.5).
 * The Python-side NativeOSXInstaller class handles Tk 8.5 by using
 * native macOS osascript dialogs instead of tkinter widgets. */
static int python_has_tkinter(const char *python_path) {
    if (access(python_path, X_OK) != 0)
        return 0;

    char cmd[1024];
    snprintf(cmd, sizeof(cmd),
             "%s -c 'import tkinter' >/dev/null 2>&1", python_path);

    int rc = system(cmd);
    return (rc == 0);
}

/* Weakest check: Python exists and is executable (no tkinter needed).
 * The Python-side code will fall back to CLI mode when tkinter is
 * unavailable, so we can still run the installer. */
static int python_is_usable(const char *python_path) {
    if (access(python_path, X_OK) != 0)
        return 0;

    char cmd[1024];
    snprintf(cmd, sizeof(cmd),
             "%s -c 'import sys; exit(0)' >/dev/null 2>&1", python_path);

    int rc = system(cmd);
    return (rc == 0);
}

/* Show a native macOS alert when no suitable Python is found. */
static void show_alert(void) {
    system(
        "osascript -e '"
        "display alert \"Python Not Found\" "
        "message \"GLINT Installer requires Python 3 to run.\\n\\n"
        "No usable Python interpreter was found on this system.\\n\\n"
        "Please install Python from python.org (3.12 / 3.13 recommended) "
        "or via Homebrew / Miniforge, then relaunch the installer.\\n\\n"
        "Details in ~/Library/Logs/GLINT_Installer.log\"'"
    );
}

/* ── Main ────────────────────────────────────────────────────────── */

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    log_open();
    log_msg("GLINT Installer launcher started");

    /* ── Set up a rich PATH so that Python discovery works
     *    even when Finder launches the .app with a minimal PATH. ── */
    const char *home = getenv("HOME");
    if (!home) home = "/tmp";

    char rich_path[4096];
    snprintf(rich_path, sizeof(rich_path),
           "/Library/Frameworks/Python.framework/Versions/Current/bin:"
           "/Library/Frameworks/Python.framework/Versions/3.13/bin:"
           "/Library/Frameworks/Python.framework/Versions/3.12/bin:"
           "/usr/local/bin:"
           "/opt/homebrew/bin:"
           "%s/miniforge3/bin:"
           "%s/mambaforge/bin:"
           "%s/miniconda3/bin:"
           "%s/anaconda3/bin:"
           "/usr/bin:/bin:/usr/sbin:/sbin",
           home, home, home, home);
    setenv("PATH", rich_path, 1);

    /* ── Find a Python interpreter with modern Tk (>= 8.6) ── */

    /* Build dynamic candidate paths for conda/miniforge/mambaforge */
    char miniforge_py[1024], mambaforge_py[1024], miniconda_py[1024], anaconda_py[1024];
    snprintf(miniforge_py,  sizeof(miniforge_py),  "%s/miniforge3/bin/python3",  home);
    snprintf(mambaforge_py, sizeof(mambaforge_py), "%s/mambaforge/bin/python3",  home);
    snprintf(miniconda_py,  sizeof(miniconda_py),  "%s/miniconda3/bin/python3",  home);
    snprintf(anaconda_py,   sizeof(anaconda_py),   "%s/anaconda3/bin/python3",   home);

    const char *candidates[] = {
        "/Library/Frameworks/Python.framework/Versions/Current/bin/python3",
        "/Library/Frameworks/Python.framework/Versions/3.13/bin/python3",
        "/Library/Frameworks/Python.framework/Versions/3.12/bin/python3",
        "/usr/local/bin/python3",
        "/opt/homebrew/bin/python3",
        miniforge_py,
        mambaforge_py,
        miniconda_py,
        anaconda_py,
        "/usr/bin/python3",
        NULL
    };

    const char *python = NULL;

    /* Pass 1: prefer a Python with Tk >= 8.6 (modern, ttk works) */
    log_msg("Pass 1: looking for Python with Tk >= 8.6 ...");
    for (int i = 0; candidates[i]; i++) {
        log_msg("Checking candidate: %s", candidates[i]);
        if (python_has_modern_tk(candidates[i])) {
            python = candidates[i];
            log_msg("  -> OK (Tk >= 8.6)");
            break;
        } else {
            log_msg("  -> skipped");
        }
    }

    /* Pass 2: fall back to any Python with tkinter (even Tk 8.5).
     * The Python-side NativeOSXInstaller class detects Tk < 8.6 and
     * uses native macOS osascript dialogs instead of tkinter widgets. */
    if (!python) {
        log_msg("Pass 2: falling back to any Python with tkinter ...");
        for (int i = 0; candidates[i]; i++) {
            log_msg("Checking candidate: %s", candidates[i]);
            if (python_has_tkinter(candidates[i])) {
                python = candidates[i];
                log_msg("  -> OK (has tkinter, Python-side will auto-select UI mode)");
                break;
            } else {
                log_msg("  -> skipped (not found or no tkinter)");
            }
        }
    }

    /* Pass 3: any usable Python at all (no tkinter required).
     * The Python-side code will fall back to CLI installer mode. */
    if (!python) {
        log_msg("Pass 3: falling back to any usable Python3 ...");
        for (int i = 0; candidates[i]; i++) {
            log_msg("Checking candidate: %s", candidates[i]);
            if (python_is_usable(candidates[i])) {
                python = candidates[i];
                log_msg("  -> OK (no tkinter, will use CLI mode)");
                break;
            } else {
                log_msg("  -> skipped (not found or not usable)");
            }
        }
    }

    if (!python) {
        log_msg("ERROR: No usable Python3 found on this system");
        show_alert();
        return 1;
    }

    log_msg("Selected Python: %s", python);

    /* ── Resolve the installer script path.
     *    This binary lives at:  <app>/Contents/MacOS/launcher
     *    The script lives at:   <app>/Contents/Resources/GLINT_Installer.py ── */

    char exe_buf[4096];
    uint32_t exe_buf_size = sizeof(exe_buf);
    if (_NSGetExecutablePath(exe_buf, &exe_buf_size) != 0) {
        log_msg("ERROR: _NSGetExecutablePath failed");
        show_alert();
        return 1;
    }

    /* Resolve symlinks so dirname() gives a canonical path. */
    char *real_exe = realpath(exe_buf, NULL);
    if (!real_exe) {
        log_msg("ERROR: realpath failed for %s", exe_buf);
        show_alert();
        return 1;
    }

    /* dirname() may modify its argument, so work on a copy. */
    char dir_copy[4096];
    strncpy(dir_copy, real_exe, sizeof(dir_copy) - 1);
    dir_copy[sizeof(dir_copy) - 1] = '\0';
    free(real_exe);

    char *macos_dir = dirname(dir_copy);          /* .../Contents/MacOS   */
    char contents_dir[4096];
    strncpy(contents_dir, macos_dir, sizeof(contents_dir) - 1);
    contents_dir[sizeof(contents_dir) - 1] = '\0';
    /* Go up one level: Contents/MacOS -> Contents */
    char *slash = strrchr(contents_dir, '/');
    if (slash) *slash = '\0';

    char installer_script[4096];
    snprintf(installer_script, sizeof(installer_script),
             "%s/Resources/GLINT_Installer.py", contents_dir);

    if (access(installer_script, R_OK) != 0) {
        log_msg("ERROR: Installer script not found: %s", installer_script);
        system(
            "osascript -e '"
            "display alert \"Installer Script Missing\" "
            "message \"GLINT_Installer.py was not found inside the app bundle.\\n"
            "The application may be corrupted. Please re-download.\"'"
        );
        return 1;
    }

    log_msg("Installer script: %s", installer_script);

    /* Log Python/Tk versions before launch. */
    {
        char ver_cmd[2048];
        snprintf(ver_cmd, sizeof(ver_cmd),
                 "%s -c '"
                 "import sys; print(\"Python:\", sys.version)\n"
                 "try:\n"
                 " import tkinter as tk; print(\"Tk:\", tk.TkVersion, \"Tcl:\", tk.TclVersion)\n"
                 "except ImportError:\n"
                 " print(\"Tk: not available\")\n"
                 "' 2>&1", python);

        FILE *fp = popen(ver_cmd, "r");
        if (fp) {
            char line[256];
            while (fgets(line, sizeof(line), fp)) {
                /* Strip trailing newline for cleaner log output. */
                size_t len = strlen(line);
                if (len > 0 && line[len - 1] == '\n') line[len - 1] = '\0';
                log_msg("  %s", line);
            }
            pclose(fp);
        }
    }

    /* Transform this process into a regular foreground application.
     * When launched from .app bundle, child processes run in Background
     * role. Tk 8.5's TkpInit crashes when calling macOS APIs that
     * require a foreground app. This must be done BEFORE launching Python. */
    ProcessSerialNumber psn = {0, kCurrentProcess};
    TransformProcessType(&psn, kProcessTransformToForegroundApplication);

    /* Remove app bundle identity and suppress Tk deprecation warning */
    unsetenv("__CFBundleIdentifier");
    setenv("TK_SILENCE_DEPRECATION", "1", 1);

    log_msg("Launching installer via execl() ...");
    if (g_log) fclose(g_log);

    /* Use execl() to completely replace this process with Python.
     * This avoids Cocoa state conflicts between the launcher and tkinter.
     * Unlike system() which forks a child that inherits the partially-
     * initialized Cocoa/NSApplication state, execl() replaces the process
     * image entirely, letting tkinter create its own NSApplication cleanly. */
    execl(python, python, installer_script, (char *)NULL);

    /* execl only returns on failure */
    perror("execl failed");
    show_alert();
    return 1;
}
