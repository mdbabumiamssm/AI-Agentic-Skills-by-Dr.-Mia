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


import asyncio

from semantic_kernel.connectors.ai.onnx import OnnxGenAIChatCompletion, OnnxGenAIPromptExecutionSettings
from semantic_kernel.contents import AuthorRole, ChatHistory, ChatMessageContent, ImageContent
from semantic_kernel.kernel import Kernel

# This concept sample shows how to use the Onnx connector with
# a local model running in Onnx

kernel = Kernel()

service_id = "phi3"
#############################################
# Make sure to download an ONNX model
# If onnxruntime-genai is used:
# (https://huggingface.co/microsoft/Phi-3-vision-128k-instruct-onnx-cpu)
# If onnxruntime-genai-cuda is installed for gpu use:
# (https://huggingface.co/microsoft/Phi-3-vision-128k-instruct-onnx-gpu)
# Then set ONNX_GEN_AI_CHAT_MODEL_FOLDER environment variable to the path to the model folder
#############################################
streaming = True

chat_completion = OnnxGenAIChatCompletion(ai_model_id=service_id, template="phi3v")

# Max length property is important to allocate RAM
# If the value is too big, you ran out of memory
# If the value is too small, your input is limited
settings = OnnxGenAIPromptExecutionSettings(max_length=4096)

system_message = """
You are a helpful assistant.
You know about provided images and the history of the conversation.
"""
chat_history = ChatHistory(system_message=system_message)


async def chat() -> bool:
    try:
        user_input = input("User:> ")
    except KeyboardInterrupt:
        print("\n\nExiting chat...")
        return False
    except EOFError:
        print("\n\nExiting chat...")
        return False

    if user_input == "exit":
        print("\n\nExiting chat...")
        return False
    chat_history.add_user_message(user_input)
    if streaming:
        print("Mosscap:> ", end="")
        message = ""
        async for chunk in chat_completion.get_streaming_chat_message_content(
            chat_history=chat_history, settings=settings, kernel=kernel
        ):
            print(chunk.content, end="")
            if chunk.content:
                message += chunk.content
        chat_history.add_assistant_message(message)
        print("")
    else:
        answer = await chat_completion.get_chat_message_content(
            chat_history=chat_history, settings=settings, kernel=kernel
        )
        print(f"Mosscap:> {answer}")
        chat_history.add_message(message)
    return True


async def main() -> None:
    chatting = True
    image_path = input("Image Path (leave empty if no image): ")
    if image_path:
        chat_history.add_message(
            ChatMessageContent(
                role=AuthorRole.USER,
                items=[
                    ImageContent.from_image_path(image_path=image_path),
                ],
            ),
        )
    while chatting:
        chatting = await chat()


if __name__ == "__main__":
    asyncio.run(main())

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
