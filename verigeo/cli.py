"""Command-Line Interface for VeriGeo (NS-VGP 1.0)."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from .protocol import (
    PROTOCOL_VERSION,
    tool_verify_and_heal_spatial_workflow,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="verigeo",
        description="VeriGeo (NS-VGP 1.0): Neuro-Symbolic Verifiable GeoAI Protocol CLI",
    )
    sub = parser.add_subparsers(dest="command")

    verify_p = sub.add_parser("verify", help="Verify and self-heal a spatial workflow step")
    verify_p.add_argument("--operator", default="spatial_join", help="Spatial operator name")
    verify_p.add_argument("--primary", default=None, help="Primary vector dataset path")
    verify_p.add_argument("--secondary", default=None, help="Secondary vector dataset path")
    verify_p.add_argument("--no-heal", action="store_true", help="Disable symbolic auto-healing")

    args = parser.parse_args(argv)
    if args.command == "verify":
        res: dict[str, Any] = tool_verify_and_heal_spatial_workflow(
            iface=None,
            args={
                "operator": args.operator,
                "primary_layer": args.primary,
                "secondary_layer": args.secondary,
                "auto_heal": not args.no_heal,
            },
        )
        print(json.dumps(res, indent=2))
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
