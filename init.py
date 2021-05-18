# -*- coding: utf-8 -*-

import src.main as main
from src.argumentParser.argumentParser import ArgumentParser

parser = ArgumentParser("Executes a parser or executes a parser and a model")

if __name__ == "__main__":
    main.run(**parser.get_function_arguments())
