import time
import os
import csv
import logging
from datetime import datetime

class Benchmark:
    def __init__(self, output_path="benchmark.tsv", enabled=False):
        self.output_path = output_path
        self.enabled = enabled
        self.timings = {}
        self.start_times = {}

    def start(self, step_name):
        """Start a timer for a specific step."""
        if not self.enabled: return
        self.start_times[step_name] = time.time()
        logging.info(f"Benchmark started for: {step_name}")

    def stop(self, step_name):
        """Stop the timer for a specific step."""
        if not self.enabled: return
        if step_name in self.start_times:
            elapsed = time.time() - self.start_times[step_name]
            # Accumulate time if step is called multiple times (though splace is linear)
            self.timings[step_name] = self.timings.get(step_name, 0) + elapsed
            logging.info(f"Benchmark finished for: {step_name} ({elapsed:.2f}s)")
            del self.start_times[step_name]
        else:
            logging.warning(f"Benchmark stop called for {step_name} without start.")

    def save(self, input_dir, file_count):
        """Save metrics to TSV."""
        if not self.enabled: return
        
        # Check if file exists to write header
        file_exists = os.path.exists(self.output_path)
        
        headers = ["Date", "Input_Directory", "File_Count", "Step", "Duration_Seconds"]
        
        try:
            with open(self.output_path, "a", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                if not file_exists:
                    writer.writerow(headers)
                
                date_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                
                for step, duration in self.timings.items():
                    writer.writerow([
                        date_str,
                        os.path.basename(os.path.abspath(input_dir)),
                        file_count,
                        step,
                        f"{duration:.4f}"
                    ])
            logging.info(f"Benchmark results saved to {self.output_path}")
        except Exception as e:
            logging.error(f"Failed to save benchmark: {e}")
