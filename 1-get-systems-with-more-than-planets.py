import pyvo
import pandas
from uio.utility.databases import tap

serviceEndpoint = tap.getServiceEndpoint("NASA")
if not serviceEndpoint:
    raise SystemError(
        "[ERROR] Couldn't get a TAP endpoint for the specified service"
    )

resultPlanets = tap.queryService(
    serviceEndpoint,
    " ".join((
        "SELECT hostname",
        "FROM ps",
        "GROUP BY hostname",
        "HAVING COUNT(DISTINCT pl_name) > 1",
        "ORDER BY hostname"
    ))
).to_table().to_pandas()

with open("./systems.txt", "w") as f:
    f.write("\n".join(resultPlanets["hostname"].values))
    f.write("\n")
