#!/usr/bin/env python3
import datetime
import os
import shutil
import stat
import subprocess
import sys


def get_target_datetime(hour=17, minute=15):
    now = datetime.datetime.now()
    target = now.replace(hour=hour, minute=minute, second=0, microsecond=0)
    if target <= now:
        target += datetime.timedelta(days=1)
    return target


def make_executable(path):
    st = os.stat(path)
    os.chmod(path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def schedule_with_at(script_path, target):
    at_bin = shutil.which("at")
    if at_bin is None:
        return False, "`at` command not found"

    when = "17:15"
    if target.date() > datetime.date.today():
        when = "17:15 tomorrow"

    command = f"exec {script_path}"
    proc = subprocess.run([at_bin, when], input=command + "\n", text=True, capture_output=True)
    if proc.returncode != 0:
        return False, proc.stderr.strip() or proc.stdout.strip()
    return True, proc.stdout.strip()


def fallback_wait_and_run(script_path, target):
    now = datetime.datetime.now()
    delay = (target - now).total_seconds()
    if delay <= 0:
        raise RuntimeError("Target time is not in the future")
    print(f"Waiting {int(delay)} seconds until {target.isoformat()}...")
    try:
        time.sleep(delay)
    except KeyboardInterrupt:
        print("Canceled.")
        sys.exit(1)
    subprocess.run([script_path], check=True)


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    target_script = os.path.join(script_dir, "run_olllp.sh")

    if not os.path.exists(target_script):
        print(f"Error: target script not found at {target_script}")
        sys.exit(1)

    make_executable(target_script)
    target_time = get_target_datetime(17, 15)

    success, message = schedule_with_at(target_script, target_time)
    if success:
        print(f"Scheduled {target_script} once at {target_time.strftime('%Y-%m-%d %H:%M')} using at.")
        if message:
            print(message)
        sys.exit(0)

    print(f"Could not schedule with at: {message}")
    print("Falling back to waiting in the current process. Keep this script running until the launch time.")
    try:
        import time
        fallback_wait_and_run(target_script, target_time)
    except Exception as exc:
        print(f"Error: {exc}")
        sys.exit(1)
