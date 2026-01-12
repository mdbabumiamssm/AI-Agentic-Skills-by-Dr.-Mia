#!/usr/bin/env python3
"""
bioskills CLI - Universal Biomedical Skills Platform
Production Version 2.0.0
"""

import sys
import os
import argparse
import yaml
import asyncio
from pathlib import Path

# Import Real Implementations
try:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from platform.biokernel.server import app, kernel
    from platform.optimizer.meta_prompter import PromptOptimizer, ModelTarget
except ImportError as e:
    print(f"Import Error: {e}")
    # Fallback for dev environment without package install
    pass

def serve_command(args):
    """Start the BioKernel Runtime Server."""
    import uvicorn
    print(f"Starting BioKernel on port {args.port}...")
    uvicorn.run("platform.biokernel.server:app", host="0.0.0.0", port=int(args.port), reload=True)
    return 0

def optimize_command(args):
    """Run the AI Prompt Optimizer."""
    input_text = args.text
    if args.file:
        with open(args.file, 'r') as f:
            input_text = f.read()
    
    if not input_text:
        print("Error: Provide text via --text or --file")
        return 1

    optimizer = PromptOptimizer()
    target = ModelTarget(args.target.lower())
    
    result = optimizer.optimize(input_text, target)
    
    print(f"\n--- Optimized for {target.name} ---\n")
    print(result)
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write(result)
        print(f"\nSaved to {args.output}")
        
    return 0

def run_agent_command(args):
    """Run a specific agent via the BioKernel."""
    # This is a CLI wrapper around the Kernel's logic
    import asyncio
    
    async def _run():
        result = await kernel.execute(args.skill, args.query, {{}})
        print(f"\nAgent Response ({result.model_used}):")
        print("-" * 40)
        print(result.response)
        print("-" * 40)
        print(f"Tools Used: {result.tools_used}")
    
    asyncio.run(_run())
    return 0

def main():
    parser = argparse.ArgumentParser(prog='bioskills')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # serve
    serve_parser = subparsers.add_parser('serve', help='Start BioKernel Server')
    serve_parser.add_argument('--port', default='8000')

    # optimize
    opt_parser = subparsers.add_parser('optimize', help='Optimize prompts')
    opt_parser.add_argument('--file', help='Input prompt file')
    opt_parser.add_argument('--text', help='Input prompt text')
    opt_parser.add_argument('--target', choices=['claude', 'openai', 'gemini'], required=True)
    opt_parser.add_argument('-o', '--output', help='Output file')

    # run-agent
    run_parser = subparsers.add_parser('run', help='Run a skill/agent')
    run_parser.add_argument('skill', help='Skill ID (e.g., prior_auth)')
    run_parser.add_argument('query', help='Input query')

    args = parser.parse_args()

    if args.command == 'serve':
        return serve_command(args)
    elif args.command == 'optimize':
        return optimize_command(args)
    elif args.command == 'run':
        return run_agent_command(args)
    else:
        parser.print_help()
        return 0

if __name__ == '__main__':
    sys.exit(main())
