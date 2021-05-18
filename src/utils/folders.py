def parse_route(route: str) -> str:
    """Ensures that a route finish with /. If not, the charaters is appended to the
    route

    Parameters
    ----------
    route: str
        Route

    Returns
    -------
    Route with / as final character
    """
    if route and not route.endswith("/"):
        route += "/"
    return route
