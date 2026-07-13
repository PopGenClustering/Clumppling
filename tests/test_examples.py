import os
from pathlib import Path
import subprocess
import sys
import zipfile

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = REPOSITORY_ROOT / "examples"


def prepare_capeverde_input(tmp_path: Path) -> Path:
    """Use an existing extracted dataset or extract the bundled ZIP."""

    extracted_input = EXAMPLES_DIR / "capeverde"

    if extracted_input.is_dir():
        return extracted_input

    archive = EXAMPLES_DIR / "capeverde.zip"
    assert archive.is_file(), f"Missing Cape Verde archive: {archive}"

    extraction_dir = tmp_path / "capeverde_extracted"
    extraction_dir.mkdir()

    with zipfile.ZipFile(archive) as zip_file:
        zip_file.extractall(extraction_dir)

    input_files = sorted(extraction_dir.rglob("*.indivq"))
    assert input_files, "No .indivq files found in capeverde.zip"

    input_directories = {path.parent for path in input_files}

    assert len(input_directories) == 1, (
        "Expected all Cape Verde .indivq files in one directory, "
        f"but found: {sorted(input_directories)}"
    )

    return input_directories.pop()


def run_capeverde(
    input_dir: Path,
    output_dir: Path,
) -> subprocess.CompletedProcess[str]:
    labels_file = EXAMPLES_DIR / "capeverde_ind_labels.txt"

    command = [
        sys.executable,
        "-m",
        "clumppling",
        "-i",
        str(input_dir),
        "-o",
        str(output_dir),
        "-f",
        "admixture",
        "--extension",
        ".indivq",
        "--ind_labels",
        str(labels_file),
        "--cd_method",
        "louvain",
        "--random_seed",
        "42",
        "-v",
        "F",
    ]

    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    environment["PYTHONHASHSEED"] = "0"

    return subprocess.run(
        command,
        cwd=REPOSITORY_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )


@pytest.mark.integration
def test_capeverde_louvain_example_is_reproducible(
    tmp_path: Path,
) -> None:
    """Run the real example twice with the same random seed."""

    input_dir = prepare_capeverde_input(tmp_path)

    output_1 = tmp_path / "capeverde_output_1"
    output_2 = tmp_path / "capeverde_output_2"

    result_1 = run_capeverde(input_dir, output_1)
    result_2 = run_capeverde(input_dir, output_2)

    assert result_1.returncode == 0, (
        f"First run failed.\n"
        f"STDOUT:\n{result_1.stdout}\n"
        f"STDERR:\n{result_1.stderr}"
    )

    assert result_2.returncode == 0, (
        f"Second run failed.\n"
        f"STDOUT:\n{result_2.stdout}\n"
        f"STDERR:\n{result_2.stderr}"
    )

    expected_directories = [
        "input",
        "alignment_withinK",
        "modes",
        "alignment_acrossK",
        "modes_aligned",
    ]

    for relative_directory in expected_directories:
        assert (output_1 / relative_directory).is_dir()
        assert (output_2 / relative_directory).is_dir()

    files_to_compare = [
        Path("alignment_acrossK/alignment_acrossK_rep.txt"),
        Path("alignment_acrossK/best_pairs_acrossK_rep.txt"),
        Path("alignment_acrossK/major_pairs_acrossK_rep.txt"),
    ]

    for relative_file in files_to_compare:
        file_1 = output_1 / relative_file
        file_2 = output_2 / relative_file

        assert file_1.is_file(), f"Missing file: {file_1}"
        assert file_2.is_file(), f"Missing file: {file_2}"

        assert file_1.read_text() == file_2.read_text(), (
            f"Seeded runs differ in {relative_file}"
        )

    mode_files_1 = sorted(
        path.relative_to(output_1)
        for path in output_1.rglob("*.Q")
    )
    mode_files_2 = sorted(
        path.relative_to(output_2)
        for path in output_2.rglob("*.Q")
    )

    assert mode_files_1
    assert mode_files_1 == mode_files_2

    for relative_file in mode_files_1:
        assert (
            output_1 / relative_file
        ).read_text() == (
            output_2 / relative_file
        ).read_text()
