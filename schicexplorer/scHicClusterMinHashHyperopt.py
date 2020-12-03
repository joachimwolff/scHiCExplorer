import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)
import argparse

from hyperopt import hp, fmin, tpe, space_eval, STATUS_OK, Trials, rand, atpe
import logging
log = logging.getLogger(__name__)
from tempfile import NamedTemporaryFile, mkdtemp

from schicexplorer import scHicClusterMinHash
from schicexplorer._version import __version__

import subprocess


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
""")

    parserRequired = parser.add_argument_group('Required arguments')

    parserRequired.add_argument('--matrix', '-m',
                                help='The matrix to compute the loops on.',
                                required=True)
    parserRequired.add_argument('--cellColor', '-c',
                                help='The file with the associated cell types or cell cycle stages.',
                                required=True)

    parserRequired.add_argument('--outputFileName', '-o',
                                help='File names for the result of the optimization'
                                ' (Default: %(default)s).',
                                default='hyperoptMinHash_result.txt',
                                required=False)
    parserOpt = parser.add_argument_group('Optional arguments')

    parserOpt.add_argument('--runs',
                           type=int,
                           default=100,
                           help='Number of runs of hyperopt.')
    parserOpt.add_argument('--nearestNeighbors', '-k',
                           type=int,
                           default=1000,
                           help='Number of runs of hyperopt.')
    parserOpt.add_argument('--numberOfHashfunctions', '-noh',
                           help='Number of hash functions range: start, stop, stepsize'
                           ' (Default: %(default)s).',
                           required=False,
                           type=int,
                           nargs=3,
                           default=(1000, 20000, 1000)
                           )
    parserOpt.add_argument('--numberOfClusters', '-noc',
                           help='Number of cluster range'
                           ' (Default: %(default)s).',
                           required=False,
                           type=int,
                           nargs=2,
                           default=(6, 15)
                           )
    parserOpt.add_argument('--numberPCADimensions', '-nop',
                           help='Number of PCA range: start, stop, stepsize'
                           ' (Default: %(default)s).',
                           required=False,
                           type=int,
                           nargs=3,
                           default=(30, 60, 1)
                           )
    parserOpt.add_argument('--umap_numberOfNeighbors', '-unon',
                           help='Number of umap neighbors range: start, stop, stepsize'
                           ' (Default: %(default)s).',
                           required=False,
                           type=int,
                           nargs=3,
                           default=(30, 60, 1)
                           )
    parserOpt.add_argument('--umap_n_components', '-unoc',
                           help='Number of umap n_components range: start, stop, stepsize'
                           ' (Default: %(default)s).',
                           required=False,
                           type=int,
                           nargs=3,
                           default=(2, 10, 1)
                           )
    parserOpt.add_argument('--umap_min_dist', '-umin',
                           help='Number of umap neighbors range: start, stop'
                           ' (Default: %(default)s).',
                           required=False,
                           type=float,
                           nargs=2,
                           default=(0.0, 0.5)
                           )
    parserRequired.add_argument('--method', '-me',
                                help='Method to optimize by hyperopt: random tree, tpe, adaptive tpe',
                                choices=['random', 'tpe', 'atpe'],
                                default='random',
                                required=True)
    parserOpt.add_argument('--threads', '-t',
                           help='Number of threads (uses the python multiprocessing module)'
                           ' (Default: %(default)s).',
                           required=False,
                           default=16,
                           type=int
                           )

    parserOpt.add_argument('--help', '-h', action='help',
                           help='Show this help message and exit.')

    parserOpt.add_argument('--version', action='version',
                           version='%(prog)s {}'.format(__version__))

    return parser


def objective(pArgs):

    if pArgs['dimensionsPCA'] <= pArgs['umap_n_components']:
        return 1

    intraChromosomalContactsOnly = ''
    noPCA = ''
    noUMAP = ''
    if pArgs['intraChromosomalContactsOnly']:
        intraChromosomalContactsOnly = '-ic'
    if pArgs['noPCA']:
        noPCA = '--noPCA'

    args = "-m {} -c {} -cm {} -o ramani_1mb_kmeans_c5.txt -t {} -k {} -nh {}  -dim_pca {} " \
        " -cct {} --colorMap glasbey_dark  {} {} {} "\
        "--umap_n_neighbors {} " \
        "--umap_n_components {}  " \
        "--umap_min_dist {}".format(
            pArgs['matrixFile'], pArgs['numberOfClusters'], pArgs['clusterMethod'], pArgs['threads'],
            pArgs['numberOfNearestNeighbors'], pArgs['numberOfHashFunctions'],
            pArgs['dimensionsPCA'], pArgs['cell_coloring_type'],
            intraChromosomalContactsOnly, noPCA, noUMAP,
            pArgs['umap_n_neighbors'], pArgs['umap_n_components'], pArgs['umap_min_dist']).split()
    scHicClusterMinHash.main(args)

    with open('correct_associated', 'r') as file:
        error_score = 1 - float(file.readline().strip())

    print('Error score: {}'.format(error_score))
    return error_score


def main(args=None):

    args = parse_arguments().parse_args(args)

    space = {

        'numberOfClusters': hp.choice('numberOfClusters', list(range(args.numberOfClusters[0], args.numberOfClusters[1], 1))),
        'clusterMethod': hp.choice('clusterMethod', ['spectral']),
        'intraChromosomalContactsOnly': hp.choice('intraChromosomalContactsOnly', [True, False]),
        'numberOfHashFunctions': hp.choice('numberOfHashFunctions', list(range(args.numberOfHashfunctions[0], args.numberOfHashfunctions[1], args.numberOfHashfunctions[2]))),
        'numberOfNearestNeighbors': args.nearestNeighbors,
        'noPCA': hp.choice('noPCA', [False, True]),
        # 'noUMAP': hp.choice('noUMAP', [True, False]),
        'dimensionsPCA': hp.choice('dimensionsPCA', list(range(args.numberPCADimensions[0], args.numberPCADimensions[1], args.numberPCADimensions[2]))),
        'umap_n_neighbors': hp.choice('umap_n_neighbors', list(range(args.umap_numberOfNeighbors[0], args.umap_numberOfNeighbors[1], args.umap_numberOfNeighbors[2]))),
        'umap_n_components': hp.choice('umap_n_components', list(range(args.umap_n_components[0], args.umap_n_components[1], args.umap_n_components[2]))),
        # 'umap_n_epochs': hp.choice('umap_n_epochs', list(range(10, 100))),
        'umap_min_dist': hp.uniform('umap_min_dist', args.umap_min_dist[0], args.umap_min_dist[1]),
        'matrixFile': args.matrix,
        'cell_coloring_type': args.cellColor,
        'threads': args.threads

    }

    # minimize the objective over the space

    trials = Trials()
    if args.method == 'random':
        best = fmin(objective, space, algo=rand.suggest, max_evals=args.runs, trials=trials)
    elif args.method == 'tpe':
        best = fmin(objective, space, algo=tpe.suggest, max_evals=args.runs, trials=trials)
    elif args.method == 'atpe':
        best = fmin(objective, space, algo=atpe.suggest, max_evals=args.runs, trials=trials)
    else:
        print('Error, hyperopt optimization method not known.')
        exit(1)

    with open(args.outputFileName, 'w') as file:
        file.write("# Created by scHiCExplorer schicClusterMinHashHyperopt {}\n\n".format(__version__))
        file.write("{}".format(space_eval(space, best)))
