# -*- coding: utf-8 -*-


from abc import ABC, abstractmethod


class AbstractModel(ABC):
    """
    This class contains default functionallity to
    generate model that will have functions to train
    and functions to get the generated model.

    This class will need data to train the model. This
    data will be tipically a list of pairs (sample, result) that
    will be used to train an specified model.

    The model will be a defined data structure, and will be created
    using a trainer.

    A trainer is a function that gets samples as input (tipically
    a list of pairs) and returns the trained model.

    The steps will be:
        - DATA -------->  PARSER --------> TRAINER --------> TRAINED_MODEL [--------> SAVER]
                                            A
                                            |
                                            |
                                           MODEL

        - DATA -------->  PARSER [--------> LOADER] --------> TRAINED_MODEL --------> RESULT
    """

    trained_model = None

    @abstractmethod
    def parser(self):
        pass

    @abstractmethod
    def trainer(self):
        pass

    @abstractmethod
    def model(self):
        pass

    @abstractmethod
    def saver(self):
        pass

    @abstractmethod
    def loader(self):
        pass
