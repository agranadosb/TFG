# -*- coding: utf-8 -*-

import src.main as main
from src.argumentParser.argumentParser import ArgumentParser
from src.model.ktssModel import KTSSModel

parser = ArgumentParser("Executes a parser or executes a parser and a model")

for i in KTSSModel.arguments:
    parser.add_argument(i)


if __name__ == "__main__":
    main.run(**parser.get_function_arguments())
