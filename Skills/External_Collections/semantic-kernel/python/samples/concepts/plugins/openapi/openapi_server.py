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

from aiohttp import web

"""OpenAPI Sample Server"""
routes = web.RouteTableDef()


@routes.post("/{name}")
async def hello(request):
    # Get path parameters
    name = request.match_info.get("name", "")
    # Get query parameters
    q = request.rel_url.query.get("q", "")
    # Get body
    body = await request.json()
    # Get headers
    headers = request.headers
    return web.Response(text=f"Hello, {name}: q={q}, body={body}, headers={headers}")


app = web.Application()
app.add_routes(routes)

if __name__ == "__main__":
    web.run_app(app)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
