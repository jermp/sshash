#!/usr/bin/env python3

import subprocess
import time
import sys
import re

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 track_rss.py <command> [args...]")
        sys.exit(1)

    cmd = sys.argv[1:]
    print(f"[TRACKER] Launching: {' '.join(cmd)}")
    
    # Start the process
    start_time = time.time()
    proc = subprocess.Popen(cmd)
    pid = proc.pid
    
    max_rss_anon_kb = 0
    
    # Regex to find the RssAnon line in the status file
    rss_anon_pattern = re.compile(r'RssAnon:\s+(\d+)\s+kB')

    # Poll the process status while it is running
    while proc.poll() is None:
        try:
            with open(f"/proc/{pid}/status", "r") as f:
                content = f.read()
                match = rss_anon_pattern.search(content)
                if match:
                    current_rss_anon = int(match.group(1))
                    if current_rss_anon > max_rss_anon_kb:
                        max_rss_anon_kb = current_rss_anon
        except FileNotFoundError:
            # The process finished and the /proc/[pid] directory is gone
            break
        except Exception as e:
            pass
        
        # Check every 0.1 seconds to catch spikes without burning CPU
        time.sleep(0.1)

    # Calculate elapsed time
    elapsed_time = time.time() - start_time
    
    # Convert peak RssAnon to Gigabytes
    max_rss_anon_gb = max_rss_anon_kb / (1024 * 1024)

    print("\n" + "="*50)
    print("   TRUE MEMORY USAGE REPORT (mmap ignored)")
    print("="*50)
    print(f"Command:        {' '.join(cmd)}")
    print(f"Wall-clock:     {elapsed_time:.2f} seconds")
    print(f"Peak RssAnon:   {max_rss_anon_kb} kB")
    print(f"Peak RssAnon:   {max_rss_anon_gb:.2f} GB")
    print("="*50 + "\n")

    sys.exit(proc.returncode)

if __name__ == "__main__":
    main()
