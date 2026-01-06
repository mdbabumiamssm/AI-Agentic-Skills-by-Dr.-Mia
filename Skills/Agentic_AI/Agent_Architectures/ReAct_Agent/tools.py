def calculator(expression):
    """Evaluates a mathematical expression."""
    try:
        # unsafe eval for demo purposes only; use a safe parser in prod
        return str(eval(expression))
    except Exception as e:
        return f"Error: {e}"

def lookup_weather(location):
    """Mock weather lookup."""
    if "london" in location.lower():
        return "Rainy, 15C"
    elif "new york" in location.lower():
        return "Sunny, 22C"
    else:
        return "Weather data unavailable for this location."

registry = {
    "Calculator": calculator,
    "Weather": lookup_weather
}
