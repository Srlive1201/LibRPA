def equal(test, reference):
    """Compare nonempty extracted text, for example required runtime flags."""
    if not test or set(test) != set(reference):
        return False, "missing extracted files"
    count = 0
    for filename in test:
        actual = test[filename]
        expected = reference[filename]
        actual = [actual] if isinstance(actual, str) else list(actual)
        expected = [expected] if isinstance(expected, str) else list(expected)
        if not actual or not expected or any(not str(x).strip() for x in actual + expected):
            return False, "empty extracted text in {}".format(filename)
        if actual != expected:
            return False, "text mismatch in {}: {} != {}".format(filename, actual, expected)
        count += len(actual)
    return True, "{} required text values match".format(count)
