# TFG


The Repository that contains the code to the mutations predictions project.


## Installation


```sh
# Clone project
git clone https://github.com/agranadosb/TFG.git

# Install poetry
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -

# Install dependences
poetry install

# If an error occurs installation of PyVCF
source .venv/bin/activate
pip3 install PyVCF
```

## Usage
```
usage: init.py [-h] [-o {p,pm}] [-p {e}] [-p_p PARSER_PREFIX] [-p_s PARSER_SUFFIX] [-m {ktss}] -vcf VCF -fasta FASTA [-s SAVE] [-r RATIO]

Executes a parser or executes a parser and a model

optional arguments:
  -h, --help            show this help message and exit
  -o {p,pm}, --operation {p,pm}
                        Operation to make: parser -> p, both -> pm
  -p {e}, --parser {e}  Parser to use: extended -> e
  -p_p PARSER_PREFIX, --parser_prefix PARSER_PREFIX
                        Length of the sequence prefix
  -p_s PARSER_SUFFIX, --parser_suffix PARSER_SUFFIX
                        Length of the sequence suffix
  -m {ktss}, --model {ktss}
                        Model to use: ktss -> KTSS_MODEL
  -vcf VCF              Route to vcf file
  -fasta FASTA          Route to fasta file
  -s SAVE, --save SAVE  Folder where save results
  -r RATIO, --ratio RATIO
                        Ratio of training and test samples
```

Ensure that before the usage you run the command:
```sh
source .venv/bin/activate
```

## Testing

```sh
source .venv/bin/activate
make test
```
