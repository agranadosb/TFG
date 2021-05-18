import random
import string


def random_string(
    size: int = 6, chars: str = string.ascii_uppercase + string.digits
) -> str:
    return "".join(random.choice(chars) for _ in range(size))
