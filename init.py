# -*- coding: utf-8 -*-

import src.main as main
import sys

if __name__ == '__main__':
    test = False
    if sys.argv[1] == '-test':
        test = True
    main.run(test)
