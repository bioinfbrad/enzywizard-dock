from __future__ import annotations
from argparse import Namespace, ArgumentParser
from ..services.dock_service import run_dock_service


def add_dock_parser(parser: ArgumentParser) -> None:
    parser.add_argument("-i", "--input_path",required=True,help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-s", "--substrate_names",required=True,help="Input substrate names separated by ','. Each substrate name must match corresponding SDF file names in substrate_dir.")
    parser.add_argument("-d","--substrate_dir",required=True,help="Path to a directory containing input substrate SDF files.")
    parser.add_argument("-o", "--output_dir",required=True,help="Directory to save docking outputs and the JSON report.")
    parser.add_argument("--max_docking_attempt_num",type=int,default=20,help="Maximum number of docking attempts (default: 20).")
    parser.add_argument("--no_early_stop", action="store_false", dest="early_stop",help="Disable stopping immediately after the first successful docking result (default: enabled).")
    parser.set_defaults(early_stop=True)
    parser.add_argument("--exhaustiveness",type=int,default=16,help="Exhaustiveness of AutoDock Vina search (default: 16). Larger values may improve docking search coverage but increase runtime.")
    parser.add_argument("--cpu",type=int,default=0,help="Number of CPUs used by AutoDock Vina (default: 0). A value of 0 lets Vina decide automatically.")
    parser.add_argument("--min_rad",type=float,default=1.8,help="Minimum probe radius used in pocket detection (default: 1.8). Smaller values may detect narrower cavities, but overly small values may cause PyVOL/MSMS failure.")
    parser.add_argument("--max_rad",type=float,default=6.2,help="Maximum probe radius used in pocket detection (default: 6.2). Larger values may detect broader cavities, but overly large values may cause PyVOL/MSMS failure.")
    parser.add_argument("--min_volume",type=int,default=50,help="Minimum pocket volume threshold used to filter detected pocket regions (default: 50). Larger values retain only larger pocket candidates.")

    parser.set_defaults(func=run_dock)


def run_dock(args: Namespace) -> None:
    run_dock_service(
        input_path=args.input_path,
        substrate_names=args.substrate_names,
        substrate_dir=args.substrate_dir,
        output_dir=args.output_dir,
        max_docking_attempt_num=args.max_docking_attempt_num,
        early_stop=args.early_stop,
        exhaustiveness=args.exhaustiveness,
        cpu=args.cpu,
        min_rad=args.min_rad,
        max_rad=args.max_rad,
        min_volume=args.min_volume,
    )


