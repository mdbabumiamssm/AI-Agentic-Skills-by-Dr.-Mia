"""
Distributed Training Simulator (2026 Update)

Uses Python's `multiprocessing` to demonstrate REAL parallel execution 
of gradient calculations, simulating a Parameter Server architecture.
"""

import time
import random
import multiprocessing
import os

def worker_task(worker_id, data_chunk, result_queue):
    """
    Simulates a training step on a GPU/Worker.
    Runs in a separate PROCESS (True Parallelism).
    """
    pid = os.getpid()
    print(f"[Worker {worker_id}] (PID {pid}) Processing {len(data_chunk)} items...")
    
    # Simulate compute heavy task
    time.sleep(1.0) 
    
    # Calculate mock gradient (average of data + noise)
    gradient = sum(data_chunk) / len(data_chunk) + random.uniform(-0.1, 0.1)
    
    print(f"[Worker {worker_id}] Finished. Gradient: {gradient:.4f}")
    result_queue.put(gradient)

class ParameterServer:
    def __init__(self):
        self.weights = 0.0
        self.learning_rate = 0.1

    def update(self, gradients):
        if not gradients:
            return
        avg_grad = sum(gradients) / len(gradients)
        self.weights -= self.learning_rate * avg_grad
        print(f"\n[Server] Updated Global Model Weights: {self.weights:.4f}")

def run_simulation():
    print(f"--- Distributed System Demo (Main PID: {os.getpid()}) ---")
    
    # Data Sharding
    dataset = list(range(100))
    # Split into 4 chunks
    num_workers = 4
    chunk_size = len(dataset) // num_workers
    batches = [dataset[i:i + chunk_size] for i in range(0, len(dataset), chunk_size)]
    
    server = ParameterServer()
    
    # Epoch Loop
    for epoch in range(1, 3):
        print(f"\n=== Epoch {epoch} ===")
        
        # Queue for inter-process communication
        result_queue = multiprocessing.Queue()
        processes = []
        
        # Scatter (Launch Processes)
        for i in range(num_workers):
            p = multiprocessing.Process(
                target=worker_task, 
                args=(i, batches[i], result_queue)
            )
            processes.append(p)
            p.start()
        
        # Gather (Wait for Processes)
        for p in processes:
            p.join()
            
        # Collect Results
        gradients = []
        while not result_queue.empty():
            gradients.append(result_queue.get())
            
        # Update Model
        server.update(gradients)

if __name__ == "__main__":
    run_simulation()