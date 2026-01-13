"""
USDL (Universal Skill Description Layer) Transpiler
---------------------------------------------------
Transforms a canonical skill specification into provider-specific
payloads so a single skill can run across OpenAI, Anthropic, Gemini, etc.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from enum import Enum
from typing import Any, Dict, List


@dataclass
class FieldSpec:
    name: str
    description: str
    type: str


@dataclass
class USDLSpec:
    """
    Canonical representation of a skill.
    """

    name: str
    description: str
    inputs: List[FieldSpec]
    outputs: List[FieldSpec]
    safety_checks: List[str]
    audit_policy: str

    @classmethod
    def from_dict(cls, payload: Dict[str, Any]) -> "USDLSpec":
        return cls(
            name=payload["name"],
            description=payload["description"],
            inputs=[FieldSpec(**field) for field in payload.get("inputs", [])],
            outputs=[FieldSpec(**field) for field in payload.get("outputs", [])],
            safety_checks=payload.get("safety_checks", []),
            audit_policy=payload.get("audit_policy", "Log decisions to BioKernel event bus."),
        )


class Provider(Enum):
    OPENAI = "openai"
    ANTHROPIC = "anthropic"
    GEMINI = "gemini"


class USDLTranspiler:
    """
    Converts USDL specs into provider-specific prompt + schema bundles.
    """

    def compile(self, spec: USDLSpec, provider: Provider) -> Dict[str, Any]:
        if provider == Provider.OPENAI:
            return self._compile_openai(spec)
        if provider == Provider.ANTHROPIC:
            return self._compile_anthropic(spec)
        if provider == Provider.GEMINI:
            return self._compile_gemini(spec)
        raise ValueError(f"Unsupported provider {provider}")

    # ------------------------------------------------------------------
    # Provider implementations
    # ------------------------------------------------------------------

    def _compile_openai(self, spec: USDLSpec) -> Dict[str, Any]:
        system = (
            f"You are the '{spec.name}' healthcare agent. "
            f"{spec.description} Always emit JSON that matches the declared schema "
            "and cite relevant clinical policies."
        )
        instructions = {
            "provider": Provider.OPENAI.value,
            "system": system,
            "input_schema": self._json_schema(spec.inputs),
            "output_schema": self._json_schema(spec.outputs),
            "safety": spec.safety_checks,
            "audit_policy": spec.audit_policy,
        }
        return instructions

    def _compile_anthropic(self, spec: USDLSpec) -> Dict[str, Any]:
        thinking_block = self._render_thinking_block(spec)
        payload = {
            "provider": Provider.ANTHROPIC.value,
            "prompt": thinking_block,
            "input_tags": self._xml_schema(spec.inputs),
            "output_tags": self._xml_schema(spec.outputs),
            "safety": spec.safety_checks,
            "audit_policy": spec.audit_policy,
        }
        return payload

    def _compile_gemini(self, spec: USDLSpec) -> Dict[str, Any]:
        prompt = (
            f"Role: {spec.name}\n"
            f"Objective: {spec.description}\n"
            "Inputs:\n"
            + "\n".join(f"- {field.name}: {field.description}" for field in spec.inputs)
            + "\nOutputs:\n"
            + "\n".join(f"- {field.name}: {field.description}" for field in spec.outputs)
            + "\nSafety:\n"
            + "\n".join(f"- {item}" for item in spec.safety_checks)
        )
        return {
            "provider": Provider.GEMINI.value,
            "prompt": prompt,
            "audit_policy": spec.audit_policy,
            "output_constraints": [field.name for field in spec.outputs],
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _json_schema(self, fields: List[FieldSpec]) -> Dict[str, Any]:
        properties = {field.name: {"type": field.type, "description": field.description} for field in fields}
        required = [field.name for field in fields]
        return {"type": "object", "properties": properties, "required": required}

    def _xml_schema(self, fields: List[FieldSpec]) -> List[str]:
        return [f"<{field.name}>{field.description}</{field.name}>" for field in fields]

    def _render_thinking_block(self, spec: USDLSpec) -> str:
        input_lines = "\n".join(f"<input name='{f.name}' type='{f.type}' />" for f in spec.inputs)
        output_lines = "\n".join(f"<output name='{f.name}' type='{f.type}' />" for f in spec.outputs)
        safety_lines = "\n".join(f"<policy>{policy}</policy>" for policy in spec.safety_checks)
        return (
            "<system>\n"
            f"<role>{spec.name}</role>\n"
            f"<description>{spec.description}</description>\n"
            "<thinking>Please reason inside <thinking> tags before responding.</thinking>\n"
            "<inputs>\n"
            f"{input_lines}\n"
            "</inputs>\n"
            "<outputs>\n"
            f"{output_lines}\n"
            "</outputs>\n"
            "<safety>\n"
            f"{safety_lines}\n"
            "</safety>\n"
            "</system>"
        )


def _demo() -> None:
    import argparse
    import json
    import sys

    parser = argparse.ArgumentParser(description="USDL Transpiler CLI")
    parser.add_argument("--file", type=str, help="Path to a USDL JSON spec file")
    parser.add_argument("--provider", type=str, choices=[p.value for p in Provider], default="openai")
    
    # If no args, run the hardcoded demo
    if len(sys.argv) == 1:
        print("Running internal demo... (Use --file <path> to transpile your own specs)")
        spec = USDLSpec(
            name="Prior Authorization Coworker",
            description="Review requests against payer policy and emit determination JSON.",
            inputs=[FieldSpec(name="clinical_note", description="Unstructured note text", type="string")],
            outputs=[FieldSpec(name="decision", description="APPROVED or DENIED", type="string")],
            safety_checks=["Never hallucinate policy references.", "Escalate unclear cases to humans."],
            audit_policy="Post determination + trace to Event Bus topic prior_auth.decisions.",
        )
        transpiler = USDLTranspiler()
        for provider in Provider:
            artifact = transpiler.compile(spec, provider)
            print(f"\n--- {provider.value.upper()} ---")
            print(artifact)
        return

    args = parser.parse_args()
    
    if args.file:
        try:
            with open(args.file, "r") as f:
                data = json.load(f)
            spec = USDLSpec.from_dict(data)
            transpiler = USDLTranspiler()
            provider_enum = Provider(args.provider)
            artifact = transpiler.compile(spec, provider_enum)
            print(json.dumps(artifact, indent=2))
        except Exception as e:
            print(f"Error processing file: {e}")
            sys.exit(1)

if __name__ == "__main__":
    _demo()
