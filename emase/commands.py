# standard library imports
from pathlib import Path
from typing import Annotated
import logging

# 3rd party library imports
from rich.logging import RichHandler
import typer

# local library imports
from emase import emase_utils


ch1 = RichHandler(level=logging.NOTSET, show_level=True, show_time=True, show_path=False, omit_repeated_times=False)
ch1.addFilter(emase_utils.NoDebugLogFilter())
ch2 = RichHandler(level=logging.NOTSET, show_level=True, show_time=True, show_path=True, omit_repeated_times=False)
ch2.addFilter(emase_utils.DebugLogFilter())

logging.basicConfig(
    level="NOTSET",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[ch1, ch2]
)

app = typer.Typer()


@app.command(help="convert a BAM file to EMASE format")
def bam2emase(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help="bam file to convert")],
    haplotypes: Annotated[list[str], typer.Option('-h', '--haplotype-char', help='haplotype, either one per -h option, i.e. -h A -h B -h C, or a shortcut -h A,B,C')],
    locusid_file: Annotated[Path, typer.Option('-m', '--locus-ids', exists=True, dir_okay=False, resolve_path=True, help='filename for the locus (usually transcripts) info')],
    out_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help="EMASE file (hdf5 format)")] = None,
    delim: Annotated[str, typer.Option('-d', '--delim', help='delimiter string between locus and haplotype in BAM file')] = '_',
    index_dtype: Annotated[str, typer.Option('--index-dtype', help='advanced option, see internal code')] = 'uint32',
    data_dtype: Annotated[str, typer.Option('--data-dtype', help='advanced_option, see internal code')] = 'uint8',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('bam2emase')
    try:
        # haplotype shortcut: the following command line options are all equal
        # -h A,B,C,D,E,F,G,H
        # -h A -h B -h C -h D -h E -h F -h G -h H
        # -h A,B,C,D -h E -h F -h G,H
        all_haplotypes: list[str] = []
        for x in haplotypes:
            all_haplotypes.extend(x.split(','))

        emase_utils.bam2emase(
            alignment_file=alignment_file,
            haplotypes=all_haplotypes,
            locusid_file=locusid_file,
            out_file=out_file,
            delim=delim,
            index_dtype=index_dtype,
            data_dtype=data_dtype
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="combine EMASE files")
def combine(
    emase_files: Annotated[list[Path], typer.Option('-i', '--emase-file', exists=True, dir_okay=False, resolve_path=True, help='EMASE file to compress, can seperate files by "," or have multiple -i')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help='name of the compressed EMASE file')],
    comp_lib: Annotated[str, typer.Option('-c', '--comp-lib', help='compression library to use')] = 'zlib',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('combine')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_emase_files: list[str] = []
        for x in emase_files:
            all_emase_files.extend(str(x).split(','))

        emase_utils.combine(
            emase_files=all_emase_files,
            output_file=output_file,
            comp_lib=comp_lib
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="count alignments")
def count_alignments(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help="alignment incidence file (h5)")],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help="gene ID to isoform ID mapping info (tsv)")],
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='basename of all the generated output files')] = 'emase',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('count_alignments')
    try:
        emase_utils.count_alignments(
            alignment_file=alignment_file,
            group_file=group_file,
            outbase=outbase
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="count shared multiread pairwise alignments")
def count_alignments(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help="EMASE file")],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help="gene ID to isoform ID mapping info (tsv)")],
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='basename of all the generated output files')] = 'emase',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('count_shared_multireads_pairwise')
    try:
        emase_utils.count_shared_multireads_pairwise(
            alignment_file=alignment_file,
            group_file=group_file,
            outbase=outbase
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="hybridize Fasta files")
def create_hybrid(
    fasta_list: Annotated[list[Path|str], typer.Option('-F', '--target-files', exists=True, dir_okay=False, resolve_path=True, help='Fasta file to parse, can seperate files by "," or have multiple -i')],
    haplotypes: Annotated[list[str], typer.Option('-s', '--suffices', help='haplotype, either one per -h option, i.e. -h A -h B -h C, or a shortcut -h A,B,C')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True)] = 'gbrs.hybridized.targets.fa',
    build_bowtie_index: bool = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('create_hybrid')
    try:
        # haplotype shortcut: the following command line options are all equal
        # -h A,B,C,D,E,F,G,H
        # -h A -h B -h C -h D -h E -h F -h G -h H
        # -h A,B,C,D -h E -h F -h G,H
        all_haplotypes: list[str] = []
        for x in haplotypes:
            all_haplotypes.extend(x.split(','))

        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_fasta_files: list[str] = []
        for x in fasta_list:
            all_fasta_files.extend(str(x).split(','))

        emase_utils.create_hybrid(
            fasta_list=all_fasta_files,
            haplotypes=all_haplotypes,
            output_file=str(output_file),
            build_bowtie_index=build_bowtie_index
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="get the common alignments")
def get_common_alignments(
    emase_files: Annotated[list[Path], typer.Option('-i', '--emase-file', exists=True, dir_okay=False, resolve_path=True, help='EMASE file to compress, can seperate files by "," or have multiple -i')],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help="EMASE file with unique reads")] = None,
    comp_lib: Annotated[str, typer.Option('-c', '--comp-lib', help='compression library to use')] = 'zlib',
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('get_common_alignments')
    try:
        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_emase_files: list[str] = []
        for x in emase_files:
            all_emase_files.extend(str(x).split(','))

        emase_utils.get_common_alignments(
            emase_files=all_emase_files,
            output_file=output_file,
            comp_lib=comp_lib
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="")
def pull_out_unique_reads(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help="Alignment profile (pseudo-alignment) in EMASE format")],
    output_file: Annotated[Path, typer.Option('-o', '--output', exists=False, dir_okay=False, writable=True, resolve_path=True, help="EMASE file with unique reads")],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help="gene ID to isoform ID mapping info (tsv)")] = None,
    shallow: Annotated[bool, typer.Option('-s', '--shallow', help='return shallow EMASE file')] = False,
    ignore_alleles: Annotated[bool, typer.Option('-a', '--ignore-alleles', help='do not require allele level uniqueness')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('pull_out_unique_reads')
    try:
        emase_utils.pull_out_unique_reads(
            alignment_file=alignment_file,
            group_file=group_file,
            output_file=output_file,
            shallow=shallow,
            ignore_alleles=ignore_alleles
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="prepare EMASE")
def prepare(
    genome_file: Annotated[list[Path], typer.Option('-G', '--genome-file', exists=True, dir_okay=False, resolve_path=True, help='Genome files, can seperate files by "," or have multiple -G')],
    haplotypes: Annotated[list[str], typer.Option('-s', '--haplotype-char', help='haplotype, either one per -h option, i.e. -h A -h B -h C, or a shortcut -h A,B,C')] = None,
    gtf_file: Annotated[list[Path], typer.Option('-g', '--gtf-file', exists=True, dir_okay=False, resolve_path=True, help='Gene Annotation File files, can seperate files by "," or have multiple -G')] = None,
    out_dir: Annotated[str, typer.Option('-o', '--out_dir', help='Output folder to store results (default: the current working directory)')] = None,
    save_g2tmap: Annotated[bool, typer.Option('-m', '--save-g2tmap', help='saves gene id to transcript id list in a tab-delimited text file')] = False,
    save_dbs: Annotated[bool, typer.Option('-d', '--save-dbs', help='save dbs')] = False,
    no_bowtie_index: Annotated[bool, typer.Option('-m', '--no-bowtie-index', help='skips building bowtie index')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('run')
    try:
        # haplotype shortcut: the following command line options are all equal
        # -h A,B,C,D,E,F,G,H
        # -h A -h B -h C -h D -h E -h F -h G -h H
        # -h A,B,C,D -h E -h F -h G,H
        all_haplotypes: list[str] = []
        for x in haplotypes:
            all_haplotypes.extend(x.split(','))

        # file shortcut: the following command line options are all equal
        # -i abc.h5 -i def.h5
        # -i abc.h5,def.h5
        all_genome_files: list[str] = []
        for x in genome_file:
            all_genome_files.extend(str(x).split(','))

        all_gtf_files: list[str] = []
        for x in gtf_file:
            all_gtf_files.extend(str(x).split(','))

        emase_utils.prepare(
            genome_files=all_genome_files,
            haplotypes=all_haplotypes,
            gtf_files=all_gtf_files,
            out_dir=out_dir,
            save_g2tmap=save_g2tmap,
            save_dbs=save_dbs,
            no_bowtie_index=no_bowtie_index
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


@app.command(help="run EMASE")
def run(
    alignment_file: Annotated[Path, typer.Option('-i', '--alignment-file', exists=True, dir_okay=False, resolve_path=True, help='EMASE alignment incidence file (in hdf5 format)')],
    group_file: Annotated[Path, typer.Option('-g', '--group-file', exists=True, dir_okay=False, resolve_path=True, help='tab delimited file of gene to transcript mapping')] = None,
    length_file: Annotated[Path, typer.Option('-L', '--length-file', exists=True, dir_okay=False, resolve_path=True, help='tab delimited file of locus(transcript) and length')] = None,
    outbase: Annotated[str, typer.Option('-o', '--outbase', help='basename of all the generated output files')] = 'emase',
    multiread_model: Annotated[int, typer.Option('-M', '--multiread-model', help='emase model (default: 4)')] = 4,
    pseudocount: Annotated[float, typer.Option('-p', '--pseudocount', help='prior read count (default: 0.0)')] = 0.0,
    read_length: Annotated[int, typer.Option('-l', '--read-length', help='specify read length')] = 100,
    max_iters: Annotated[int, typer.Option('-m', '--max-iters', help='maximum iterations for EM iteration')] = 999,
    tolerance: Annotated[float, typer.Option('-t', '--tolerance', help='tolerance for EM termination (default: 0.0001 in TPM)')] = 0.0001,
    report_alignment_counts: Annotated[bool, typer.Option('-c', '--report-alignment-counts', help='whether to report alignment counts')] = False,
    report_posterior: Annotated[bool, typer.Option('-w', '--report-alignment-counts', help='whether to report posterior probabilities')] = False,
    verbose: Annotated[int, typer.Option('-v', '--verbose', count=True, help="specify multiple times for more verbose output")] = 0
) -> None:
    logger = emase_utils.configure_logging(verbose)
    logger.debug('run')
    try:
        if not multiread_model in (1, 2, 3, 4):
            raise typer.Abort(f'-M, --multiread-model must be one of 1, 2, 3, or 4')

        emase_utils.run(
            alignment_file=alignment_file,
            group_file=group_file,
            length_file=length_file,
            outbase=outbase,
            multiread_model=multiread_model,
            read_length=read_length,
            pseudocount=pseudocount,
            max_iters=max_iters,
            tolerance=tolerance,
            report_alignment_counts=report_alignment_counts,
            report_posterior=report_posterior
        )
    except Exception as e:
        if logger.level == logging.DEBUG:
            logger.exception(e)
        else:
            logger.error(e)


if __name__ == '__main__':
    app()
