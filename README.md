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
usage: init.py [-h] [-m {ktss}] [-o {p,pm}] [-p {e,m}] -s SAVE [-p_p PARSER_PREFIX] [-p_s PARSER_SUFFIX] [-r RATIO] -vcf VCF -fasta FASTA [-k K] [-ktss_nas] [-amto] [-ao] [-wc] [-pfilename PARSER_FILENAME] [-sep SEPARATOR] [-min] [-aoval] [-amv]

Executes a parser or executes a parser and a model

optional arguments:
  -h, --help            show this help message and exit
  -m {ktss}, --model {ktss}
                        Model to use: ktss -> ktss
  -o {p,pm}, --operation {p,pm}
                        Operation to make: parser -> p, both -> pm
  -p {e,m}, --parser {e,m}
                        Parser to use: extended -> e, mutation type -> m
  -s SAVE, --save SAVE  Folder where save results
  -p_p PARSER_PREFIX, --parser_prefix PARSER_PREFIX
                        Length of the sequence prefix
  -p_s PARSER_SUFFIX, --parser_suffix PARSER_SUFFIX
                        Length of the sequence suffix
  -r RATIO, --ratio RATIO
                        Ratio of training and test samples
  -vcf VCF, --vcf VCF   
                        Route to vcf file
  -fasta FASTA, --fasta FASTA
                        Route to fasta file
  -k K, --k K           
                        k value for ktss model
  -ktss_nas, --ktss-not-allowed-segments
                        Create not allowed segments
  -amto, --add-mutation-to-original
                        Add mutation to original sequence on parser file
  -ao, --add-original   
                        Writes original sequence into parsed sequences file on parser file
  -wc, --write-chromosme
                        Writes the chromsome where the sequence are from on parser file
  -pfilename PARSER_FILENAME, --parser-filename PARSER_FILENAME
                        Filename of the parser file
  -sep SEPARATOR, --separator SEPARATOR
                        Specifies the separator between characters of each sequence of the valdiator
  -min, --minimum       
                        If true only returns the minimum value and infix of the all the distances of the valdiator
  -aoval, --add-original-validator
                        If true returns the original sequence, if not returns the anotated sequence of the valdiator
  -amv, --add-mutation-validator
                        If true returns the mutation, if not returns the reference sequence
```

Ensure that before the usage you run the command:
```sh
source .venv/bin/activate
```

## Docs
```sh
# To generate the docs
make build-docs
```

## Testing

```sh
source .venv/bin/activate
make test
```
