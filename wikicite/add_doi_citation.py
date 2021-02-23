from wikidataintegrator import wdi_core, wdi_login, wdi_helpers
import os

if "WDUSER" in os.environ and "WDPASS" in os.environ and "DOI" in os.environ:
    WDUSER = os.environ['WDUSER']
    WDPASS = os.environ['WDPASS']
    DOI = os.environ["DOI"]
else:
    raise ValueError("WDUSER and WDPASS must be specified in  as environment variables. e.g. docker run -it -e WDUSER=<Wikidata username> -e WDPASS=<Wikidata password> -e DOI=10.1186/s12915-020-00940-y micelio/wikicite_align_citation:latest")

login = wdi_login.WDLogin(WDUSER, WDPASS)

qid = wdi_helpers.PublicationHelper(DOI, id_type="doi",
                                                     source="crossref").get_or_create(login if True else None)

print(qid)