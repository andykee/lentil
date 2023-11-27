# ref: https://gist.github.com/peterhurford/09f7dcda0ab04b95c026c60fa49c2a68?permalink_comment_id=3453153#gistcomment-3453153

from glob import glob


def refactor(string: str) -> str:
    return string.replace("/", ".").replace("\\", ".").replace(".py", "")


pytest_plugins = [
    refactor(fixture) for fixture in glob("tests/fixtures/*.py") if "__" not in fixture
]
