import argparse
import pprint
import re
import sys
from matplotlib import pyplot
import seaborn
import random
from collections import defaultdict
from matplotlib.ticker import StrMethodFormatter, MaxNLocator
from RW_group import Group
from RW_strengths import Strengths, History
from itertools import zip_longest

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Behold! My Rescorla-Wagnerinator!",
        epilog = '--alpha_[A-Z] ALPHA\tAssociative strength of CS A..Z. By default 0',
    )

    parser.add_argument("--beta", type = float, default = 1, help="Associativity of the US +.")
    parser.add_argument("--beta-neg", type = float, default = None, help="Associativity of the absence of US +. Equal to beta by default.")
    parser.add_argument("--lamda", type = float, default = 1, help="Asymptote of learning.")

    parser.add_argument("--use-configurals", type = bool, action = argparse.BooleanOptionalAction, help = 'Use compound stimuli with configural cues')

    parser.add_argument("--adaptive-type", choices = ['linear', 'exponential', 'macknhall'], help = 'Type of adaptive attention mode to use')
    parser.add_argument("--window-size", type = int, default = None, help = 'Size of sliding window for adaptive learning')

    parser.add_argument("--xi-hall", type = float, default = 0.2, help = 'Xi parameter for Hall alpha calculation')

    parser.add_argument("--plot-experiments", nargs = '*', help = 'List of experiments to plot. By default plot everything')
    parser.add_argument("--plot-stimuli", nargs = '*', help = 'List of stimuli, compound and simple, to plot. By default plot everything')
    parser.add_argument('--plot-alphas', type = bool, action = argparse.BooleanOptionalAction, help = 'Whether to plot the alphas of all CS.')

    parser.add_argument(
        "experiment_file",
        nargs='?',
        type = argparse.FileType('r'),
        default = sys.stdin,
        help="Path to the experiment file."
    )

    # Also accept arguments of the form --alpha_[A-Z]=n.
    args, rest = parser.parse_known_args()
    args.alphas = dict()
    for arg in rest:
        match = re.fullmatch('--alpha_([A-Z])\s*=?\s*([0-9]*\.?[0-9]*)', arg)
        if not match:
            parser.error(f'Option not understood: {arg}')

        args.alphas[match.group(1)] = float(match.group(2))

    if args.beta_neg is None:
        args.beta_neg = args.beta

    if args.use_configurals is None:
        args.use_configurals = False

    args.use_adaptive = args.adaptive_type is not None

    if args.adaptive_type == 'macknhall' and args.window_size is None:
        args.window_size = 3

    return args

def parse_parts(phase : str) -> tuple[list[tuple[str, str]], None | float]:
    rand = False
    lamda = None

    parts = []
    for part in phase.strip().split('/'):
        if part == 'rand':
            rand = True
        elif (match := re.fullmatch(r'lamb?da *= *([0-9]*(?:\.[0-9]*)?)', part)) is not None:
            lamda = float(match.group(1))
        elif (match := re.fullmatch(r'([0-9]*)([A-Z]+)([+-]?)', part)) is not None:
            num, cs, sign = match.groups()
            parts += int(num or '1') * [(cs, sign or '+')]
        else:
            raise ValueError(f'Part not understood: {part}')

    if rand:
        random.shuffle(parts)

    return parts, lamda

def run_group_experiments(g, experiment : list[tuple[list[tuple[str, str]], None | float]]) -> list[list[Strengths]]:
    results = []

    for trial, (phase, phase_lamda) in enumerate(experiment):
        strength_hist = g.runPhase(phase, phase_lamda)
        results.append(strength_hist)

    return results

def plot_graphs(data: list[dict[str, History]], plot_alphas = False):
    seaborn.set()
    pyplot.ion()

    for phase_num, experiments in enumerate(data):
        if not plot_alphas:
            fig, axes = pyplot.subplots(1, 1, figsize = (8, 6))
            axes = [axes]
        else:
            fig, axes = pyplot.subplots(1, 2, figsize = (16, 6))

        colors = dict(zip(experiments.keys(), seaborn.color_palette('husl', len(experiments))))
        for key, hist in experiments.items():
            axes[0].plot(hist.assoc, label=key, marker='D', color = colors[key], markersize=4, alpha=.5)

            if plot_alphas:
                axes[1].plot(hist.alpha, label=key + r' - $\alpha$', color = colors[key], marker='D', markersize=8, alpha=.5)
                axes[1].plot(hist.alpha_mack, label=key + r' - $\alpha_{MACK}$', color = colors[key], marker='$M$', markersize=8, alpha=.5)
                axes[1].plot(hist.alpha_hall, label=key + r' - $\alpha_{HALL}$', color = colors[key], marker='$H$', markersize=8, alpha=.5)

        axes[0].set_xlabel('Trial Number')
        axes[0].set_ylabel('Associative Strength')
        axes[0].set_title(f'Phase {phase_num} Associative Strengths')
        axes[0].legend()

        if plot_alphas:
            axes[1].set_xlabel('Trial Number')
            axes[1].set_ylabel('Alpha')
            axes[1].set_title(f'Phase {phase_num} Alphas')
            axes[1].legend()

        pyplot.tight_layout()
        pyplot.show()

    input('Press any key to continue...')


def main():
    args = parse_args()

    groups_strengths = []

    for e, experiment in enumerate(args.experiment_file.readlines()):
        if experiment not in (args.plot_experiments or [experiment]):
            continue

        name, *phases = experiment.split('|')
        name = name.strip()
        phases = [parse_parts(phase) for phase in phases]

        cs = set(''.join(y[0] for x in phases for y in x[0]))
        g = Group(name, args.alphas, args.beta_neg, args.beta, args.lamda, cs, args.use_configurals, args.adaptive_type, args.window_size, args.xi_hall)

        for phase_num, strength_hist in enumerate(run_group_experiments(g, phases)):
            while len(groups_strengths) <= phase_num:
                groups_strengths.append(defaultdict(lambda: History()))

            for strengths in strength_hist:
                for cs in strengths.combined_cs():
                    if cs not in (args.plot_stimuli or [cs]):
                        continue

                    groups_strengths[phase_num][f'{name} - {cs}'].add(strengths[cs])

    plot_graphs(groups_strengths, args.plot_alphas)

if __name__ == '__main__':
    main()
