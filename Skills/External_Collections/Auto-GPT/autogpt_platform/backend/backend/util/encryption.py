# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import json
from typing import Optional

from cryptography.fernet import Fernet

from backend.util.settings import Settings

ENCRYPTION_KEY = Settings().secrets.encryption_key


class JSONCryptor:
    def __init__(self, key: Optional[str] = None):
        # Use provided key or get from environment
        self.key = key or ENCRYPTION_KEY
        if not self.key:
            raise ValueError(
                "Encryption key must be provided or set in ENCRYPTION_KEY environment variable"
            )
        self.fernet = Fernet(
            self.key.encode() if isinstance(self.key, str) else self.key
        )

    def encrypt(self, data: dict) -> str:
        """Encrypt dictionary data to string"""
        json_str = json.dumps(data)
        encrypted = self.fernet.encrypt(json_str.encode())
        return encrypted.decode()

    def decrypt(self, encrypted_str: str) -> dict:
        """Decrypt string to dictionary"""
        if not encrypted_str:
            return {}
        try:
            decrypted = self.fernet.decrypt(encrypted_str.encode())
            return json.loads(decrypted.decode())
        except Exception:
            return {}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
