def test_imports():
    try:
        import franklab_mountainsort
        assert True
    except ImportError:
        assert False
