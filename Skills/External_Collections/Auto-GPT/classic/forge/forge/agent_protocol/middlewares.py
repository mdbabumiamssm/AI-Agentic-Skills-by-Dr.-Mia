# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from starlette.types import ASGIApp


class AgentMiddleware:
    """
    Middleware that injects the agent instance into the request scope.
    """

    def __init__(self, app: ASGIApp, agent):
        """

        Args:
            app: The FastAPI app - automatically injected by FastAPI.
            agent: The agent instance to inject into the request scope.

        Examples:
            >>> from fastapi import FastAPI, Request
            >>> from agent_protocol.agent import Agent
            >>> from agent_protocol.middlewares import AgentMiddleware
            >>> app = FastAPI()
            >>> @app.get("/")
            >>> async def root(request: Request):
            >>>     agent = request["agent"]
            >>>     task = agent.db.create_task("Do something.")
            >>>     return {"task_id": a.task_id}
            >>> agent = Agent()
            >>> app.add_middleware(AgentMiddleware, agent=agent)
        """
        self.app = app
        self.agent = agent

    async def __call__(self, scope, receive, send):
        scope["agent"] = self.agent
        await self.app(scope, receive, send)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
