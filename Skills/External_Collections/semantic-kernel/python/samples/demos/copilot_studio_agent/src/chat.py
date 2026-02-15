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

import logging
import os

import chainlit as cl
from direct_line_agent import DirectLineAgent
from dotenv import load_dotenv

from semantic_kernel.contents.chat_history import ChatHistory

load_dotenv(override=True)

logging.basicConfig(level=logging.INFO)
logging.getLogger("direct_line_agent").setLevel(logging.DEBUG)
logger = logging.getLogger(__name__)

agent = DirectLineAgent(
    id="copilot_studio",
    name="copilot_studio",
    description="copilot_studio",
    bot_secret=os.getenv("BOT_SECRET"),
    bot_endpoint=os.getenv("BOT_ENDPOINT"),
)


@cl.on_chat_start
async def on_chat_start():
    cl.user_session.set("chat_history", ChatHistory())


@cl.on_message
async def on_message(message: cl.Message):
    chat_history: ChatHistory = cl.user_session.get("chat_history")

    chat_history.add_user_message(message.content)

    response = await agent.get_response(history=chat_history)

    cl.user_session.set("chat_history", chat_history)

    logger.info(f"Response: {response}")

    await cl.Message(content=response.content, author=agent.name).send()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
