# -*- coding: utf-8 -*-

import src.main as main
from src.argumentParser.argumentParser import ArgumentParser
from src.model.ktssModel import KTSSModel
from src.model.ktssValidation import KTSSValidator
from src.parser.parserVcf import ParserVcf

parser = ArgumentParser("Executes a parser or executes a parser and a model")

""" TODO: Hacer una clase Runner que haga todo esto instanciandola """
""" TODO: Hacer que se detecten automáticamente todos las clases de argumentos """
for i in KTSSModel.arguments + ParserVcf.arguments + KTSSValidator.arguments:
    parser.add_argument(i)


if __name__ == "__main__":
    main.run(**parser.get_function_arguments())
