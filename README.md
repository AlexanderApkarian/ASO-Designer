# ASO Designer - by Alexander Apkarian

Standalone Python version of the ASO Design Template workbook.

It keeps the spreadsheet's input model:

- three flexible identifiers, for example gene, variant, and chemical modification pattern
- ASO chemistry preset or custom per-position chemistry pattern
- variant type, start position (first base = 1), length, and step size
- RNA sequence entered 5' to 3'; spaces, hyphens, and underscores can be used as visual gap/deletion markers and are ignored

It produces the same core outputs:

- ASO microwalk, advancing one base per row
- IDT ordering code
- display sequence with insertion/substitution bases highlighted or deletion gaps shown
- alignment grid equivalent to the spreadsheet output from row 42 downward
- Excel `.xlsx` export

## Install

Use Python 3.10 or newer.

```bash
python3 -m pip install -r requirements.txt
```

Tkinter is included with most Python desktop installs. If your Python build does not include Tkinter, install a standard Python.org build on macOS or enable your OS package for Tkinter.

## Run The Desktop App

```bash
python3 aso_designer.py
```

Enter the sequence and inputs, click **Calculate**, then **Export Excel**.

## Command-Line Export

```bash
python3 aso_designer.py \
  --export aso_output.xlsx \
  --target-gene GENE \
  --snp-identifier Example \
  --chemistry-number C1 \
  --mutation-type Insertion \
  --mutation-length 1 \
  --mutation-start 22 \
  --rna-sequence AUGCUACGUAUGCUACGUAUGGCAUCGUAUGCUACGUAUGCUACGUA
```

For long sequences, put the RNA in a text file and use:

```bash
python3 aso_designer.py --export aso_output.xlsx --rna-file sequence.txt --mutation-type Deletion --mutation-length 4 --mutation-start 22
```

## Notes

- The original workbook reverses the clean RNA sequence before microwalk generation. This app mirrors that behavior.
- If the supplied RNA is too short to generate the complete walk around the variant, the app returns every ASO it can design and flags the result as a partial walk.
- IDT code generation supports the workbook gapmer presets plus per-position custom base modifications/linkages.
- Choose `Custom` to open the bubble editor for DNA, DNA + 5MeC, LNA, MOE, 2'OMe, 2'F, and PS/PO linkages.
- The KT777/valeriasen preset uses a 5-10-5 MOE/DNA pattern with mixed PS/PO linkages and 5MeC C bases.
