"""ppi-hub-finder: protein interaction lookups (STRING, BioGRID) and hub-target ranking."""
import warnings

# macOS system Python links LibreSSL, and urllib3 v2 warns about it on import. It is about
# the interpreter's TLS library, not this package, and every worker process would repeat it.
warnings.filterwarnings("ignore", message=r".*OpenSSL.*", module="urllib3.*")
