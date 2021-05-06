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

## Testing

```sh
source .venv/bin/activate
make test
```
