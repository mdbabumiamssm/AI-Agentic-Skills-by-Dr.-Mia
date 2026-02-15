# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

import os
import wave
from typing import ClassVar

import keyboard
import pyaudio
from pydantic import BaseModel


class AudioRecorder(BaseModel):
    """A class to record audio from the microphone and save it to a WAV file.

    To start recording, press the spacebar. To stop recording, release the spacebar.

    To use as a context manager, that automatically removes the output file after exiting the context:
    ```
    with AudioRecorder(output_filepath="output.wav") as recorder:
        recorder.start_recording()
        # Do something with the recorded audio
        ...
    ```
    """

    # Audio recording parameters
    FORMAT: ClassVar[int] = pyaudio.paInt16
    CHANNELS: ClassVar[int] = 1
    RATE: ClassVar[int] = 44100
    CHUNK: ClassVar[int] = 1024

    output_filepath: str

    def start_recording(self) -> None:
        # Wait for the spacebar to be pressed to start recording
        keyboard.wait("space")

        # Start recording
        audio = pyaudio.PyAudio()
        stream = audio.open(
            format=self.FORMAT,
            channels=self.CHANNELS,
            rate=self.RATE,
            input=True,
            frames_per_buffer=self.CHUNK,
        )
        frames = []

        while keyboard.is_pressed("space"):
            data = stream.read(self.CHUNK)
            frames.append(data)

        # Recording stopped as the spacebar is released
        stream.stop_stream()
        stream.close()

        # Save the recorded data as a WAV file
        with wave.open(self.output_filepath, "wb") as wf:
            wf.setnchannels(self.CHANNELS)
            wf.setsampwidth(audio.get_sample_size(self.FORMAT))
            wf.setframerate(self.RATE)
            wf.writeframes(b"".join(frames))

        audio.terminate()

    def remove_output_file(self) -> None:
        os.remove(self.output_filepath)

    def __enter__(self) -> "AudioRecorder":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.remove_output_file()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
