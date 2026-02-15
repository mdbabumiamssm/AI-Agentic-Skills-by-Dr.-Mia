# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""A simple spinner module"""
import itertools
import sys
import threading
import time


class Spinner:
    """A simple spinner class"""

    def __init__(
        self,
        message: str = "Loading...",
        delay: float = 0.1,
        plain_output: bool = False,
    ) -> None:
        """Initialize the spinner class

        Args:
            message (str): The message to display.
            delay (float): The delay between each spinner update.
            plain_output (bool): Whether to display the spinner or not.
        """
        self.plain_output = plain_output
        self.spinner = itertools.cycle(["-", "/", "|", "\\"])
        self.delay = delay
        self.message = message
        self.running = False
        self.spinner_thread = None

    def spin(self) -> None:
        """Spin the spinner"""
        if self.plain_output:
            self.print_message()
            return
        while self.running:
            self.print_message()
            time.sleep(self.delay)

    def print_message(self):
        sys.stdout.write(f"\r{' ' * (len(self.message) + 2)}\r")
        sys.stdout.write(f"{next(self.spinner)} {self.message}\r")
        sys.stdout.flush()

    def start(self):
        self.running = True
        self.spinner_thread = threading.Thread(target=self.spin)
        self.spinner_thread.start()

    def stop(self):
        self.running = False
        if self.spinner_thread is not None:
            self.spinner_thread.join()
        sys.stdout.write(f"\r{' ' * (len(self.message) + 2)}\r")
        sys.stdout.flush()

    def __enter__(self):
        """Start the spinner"""
        self.start()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback) -> None:
        """Stop the spinner

        Args:
            exc_type (Exception): The exception type.
            exc_value (Exception): The exception value.
            exc_traceback (Exception): The exception traceback.
        """
        self.stop()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
