import sys

try:
    from ifblup_cli import main
except ImportError:  # imported as a package submodule
    from .ifblup_cli import main

if __name__ == '__main__':
    sys.exit(main())
