# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import tweepy

from backend.blocks.twitter._auth import (
    TEST_CREDENTIALS,
    TEST_CREDENTIALS_INPUT,
    TWITTER_OAUTH_IS_CONFIGURED,
    TwitterCredentials,
    TwitterCredentialsField,
    TwitterCredentialsInput,
)
from backend.blocks.twitter.tweepy_exceptions import handle_tweepy_exception
from backend.data.block import (
    Block,
    BlockCategory,
    BlockOutput,
    BlockSchemaInput,
    BlockSchemaOutput,
)
from backend.data.model import SchemaField


class TwitterHideReplyBlock(Block):
    """
    Hides a reply of one of your tweets
    """

    class Input(BlockSchemaInput):
        credentials: TwitterCredentialsInput = TwitterCredentialsField(
            ["tweet.read", "tweet.moderate.write", "users.read", "offline.access"]
        )

        tweet_id: str = SchemaField(
            description="ID of the tweet reply to hide",
            placeholder="Enter tweet ID",
        )

    class Output(BlockSchemaOutput):
        success: bool = SchemaField(description="Whether the operation was successful")

    def __init__(self):
        super().__init__(
            id="07d58b3e-a630-11ef-a030-93701d1a465e",
            description="This block hides a reply to a tweet.",
            categories={BlockCategory.SOCIAL},
            input_schema=TwitterHideReplyBlock.Input,
            output_schema=TwitterHideReplyBlock.Output,
            disabled=not TWITTER_OAUTH_IS_CONFIGURED,
            test_input={
                "tweet_id": "1234567890",
                "credentials": TEST_CREDENTIALS_INPUT,
            },
            test_credentials=TEST_CREDENTIALS,
            test_output=[
                ("success", True),
            ],
            test_mock={"hide_reply": lambda *args, **kwargs: True},
        )

    @staticmethod
    def hide_reply(
        credentials: TwitterCredentials,
        tweet_id: str,
    ):
        try:
            client = tweepy.Client(
                bearer_token=credentials.access_token.get_secret_value()
            )

            client.hide_reply(id=tweet_id, user_auth=False)

            return True

        except tweepy.TweepyException:
            raise

    async def run(
        self,
        input_data: Input,
        *,
        credentials: TwitterCredentials,
        **kwargs,
    ) -> BlockOutput:
        try:
            success = self.hide_reply(
                credentials,
                input_data.tweet_id,
            )
            yield "success", success
        except Exception as e:
            yield "error", handle_tweepy_exception(e)


class TwitterUnhideReplyBlock(Block):
    """
    Unhides a reply to a tweet
    """

    class Input(BlockSchemaInput):
        credentials: TwitterCredentialsInput = TwitterCredentialsField(
            ["tweet.read", "tweet.moderate.write", "users.read", "offline.access"]
        )

        tweet_id: str = SchemaField(
            description="ID of the tweet reply to unhide",
            placeholder="Enter tweet ID",
        )

    class Output(BlockSchemaOutput):
        success: bool = SchemaField(description="Whether the operation was successful")

    def __init__(self):
        super().__init__(
            id="fcf9e4e4-a62f-11ef-9d85-57d3d06b616a",
            description="This block unhides a reply to a tweet.",
            categories={BlockCategory.SOCIAL},
            input_schema=TwitterUnhideReplyBlock.Input,
            output_schema=TwitterUnhideReplyBlock.Output,
            disabled=not TWITTER_OAUTH_IS_CONFIGURED,
            test_input={
                "tweet_id": "1234567890",
                "credentials": TEST_CREDENTIALS_INPUT,
            },
            test_credentials=TEST_CREDENTIALS,
            test_output=[
                ("success", True),
            ],
            test_mock={"unhide_reply": lambda *args, **kwargs: True},
        )

    @staticmethod
    def unhide_reply(
        credentials: TwitterCredentials,
        tweet_id: str,
    ):
        try:
            client = tweepy.Client(
                bearer_token=credentials.access_token.get_secret_value()
            )

            client.unhide_reply(id=tweet_id, user_auth=False)

            return True

        except tweepy.TweepyException:
            raise

    async def run(
        self,
        input_data: Input,
        *,
        credentials: TwitterCredentials,
        **kwargs,
    ) -> BlockOutput:
        try:
            success = self.unhide_reply(
                credentials,
                input_data.tweet_id,
            )
            yield "success", success
        except Exception as e:
            yield "error", handle_tweepy_exception(e)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
