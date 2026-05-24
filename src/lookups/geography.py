"""Country → continent → hemisphere lookup.

NCBI's geo_loc_name uses ISO country names but with some legacy variants
("USA" not "United States", "UK" sometimes, etc.). We normalize via a
synonym map to canonical ISO names, then look up continent + hemisphere.

This is intentionally not a complete world atlas — coverage is the
~200 countries that actually appear in NCBI Pathogen Detection. Edge
cases (territories, disputed regions, "Yugoslavia" in old records)
return None and are surfaced as "unknown" in the spread badges.

Hemisphere here means the equator-bisecting hemisphere (Northern vs.
Southern) because that's the biologically relevant axis for foodborne
disease seasonality. Countries straddling the equator (Indonesia,
Brazil, Kenya, Ecuador, ...) are classified by their population centroid.
"""

from __future__ import annotations

__all__ = [
    "canonical_country",
    "continent_of",
    "hemisphere_of",
    "is_us_state",
]


# Synonym → canonical name.
# Add aliases that actually appear in NCBI data; keep this list pruned.
_SYNONYMS = {
    "usa": "United States",
    "united states of america": "United States",
    "u.s.a.": "United States",
    "u.s.": "United States",
    "us": "United States",
    "uk": "United Kingdom",
    "great britain": "United Kingdom",
    "britain": "United Kingdom",
    "england": "United Kingdom",
    "scotland": "United Kingdom",
    "wales": "United Kingdom",
    "northern ireland": "United Kingdom",
    "viet nam": "Vietnam",
    "côte d'ivoire": "Cote d'Ivoire",
    "ivory coast": "Cote d'Ivoire",
    "south korea": "Republic of Korea",
    "korea, republic of": "Republic of Korea",
    "north korea": "North Korea",
    "russia": "Russian Federation",
    "iran": "Iran",  # NCBI inconsistent
    "syria": "Syrian Arab Republic",
    "tanzania": "Tanzania",
    "czech republic": "Czechia",
    "macedonia": "North Macedonia",
    "burma": "Myanmar",
    "east timor": "Timor-Leste",
    "cape verde": "Cabo Verde",
    "swaziland": "Eswatini",
    "the netherlands": "Netherlands",
    "holland": "Netherlands",
    "republic of ireland": "Ireland",
}


# Canonical name → (continent, hemisphere)
# Hemisphere: 'N' for northern, 'S' for southern; equator-straddlers
# assigned by population centroid.
_COUNTRY_DATA: dict[str, tuple[str, str]] = {
    # North America
    "United States":            ("North America", "N"),
    "Canada":                   ("North America", "N"),
    "Mexico":                   ("North America", "N"),
    "Cuba":                     ("North America", "N"),
    "Dominican Republic":       ("North America", "N"),
    "Haiti":                    ("North America", "N"),
    "Jamaica":                  ("North America", "N"),
    "Puerto Rico":              ("North America", "N"),
    "Trinidad and Tobago":      ("North America", "N"),
    "Bahamas":                  ("North America", "N"),
    "Barbados":                 ("North America", "N"),
    "Costa Rica":               ("North America", "N"),
    "El Salvador":              ("North America", "N"),
    "Guatemala":                ("North America", "N"),
    "Honduras":                 ("North America", "N"),
    "Nicaragua":                ("North America", "N"),
    "Panama":                   ("North America", "N"),
    "Belize":                   ("North America", "N"),

    # South America
    "Brazil":                   ("South America", "S"),
    "Argentina":                ("South America", "S"),
    "Chile":                    ("South America", "S"),
    "Peru":                     ("South America", "S"),
    "Colombia":                 ("South America", "N"),
    "Venezuela":                ("South America", "N"),
    "Ecuador":                  ("South America", "S"),
    "Bolivia":                  ("South America", "S"),
    "Paraguay":                 ("South America", "S"),
    "Uruguay":                  ("South America", "S"),
    "Guyana":                   ("South America", "N"),
    "Suriname":                 ("South America", "N"),
    "French Guiana":            ("South America", "N"),

    # Europe
    "United Kingdom":           ("Europe", "N"),
    "Germany":                  ("Europe", "N"),
    "France":                   ("Europe", "N"),
    "Italy":                    ("Europe", "N"),
    "Spain":                    ("Europe", "N"),
    "Portugal":                 ("Europe", "N"),
    "Netherlands":              ("Europe", "N"),
    "Belgium":                  ("Europe", "N"),
    "Switzerland":              ("Europe", "N"),
    "Austria":                  ("Europe", "N"),
    "Sweden":                   ("Europe", "N"),
    "Norway":                   ("Europe", "N"),
    "Denmark":                  ("Europe", "N"),
    "Finland":                  ("Europe", "N"),
    "Iceland":                  ("Europe", "N"),
    "Ireland":                  ("Europe", "N"),
    "Poland":                   ("Europe", "N"),
    "Czechia":                  ("Europe", "N"),
    "Slovakia":                 ("Europe", "N"),
    "Hungary":                  ("Europe", "N"),
    "Romania":                  ("Europe", "N"),
    "Bulgaria":                 ("Europe", "N"),
    "Greece":                   ("Europe", "N"),
    "Croatia":                  ("Europe", "N"),
    "Slovenia":                 ("Europe", "N"),
    "Serbia":                   ("Europe", "N"),
    "Bosnia and Herzegovina":   ("Europe", "N"),
    "North Macedonia":          ("Europe", "N"),
    "Albania":                  ("Europe", "N"),
    "Montenegro":               ("Europe", "N"),
    "Kosovo":                   ("Europe", "N"),
    "Ukraine":                  ("Europe", "N"),
    "Belarus":                  ("Europe", "N"),
    "Lithuania":                ("Europe", "N"),
    "Latvia":                   ("Europe", "N"),
    "Estonia":                  ("Europe", "N"),
    "Moldova":                  ("Europe", "N"),
    "Russian Federation":       ("Europe", "N"),
    "Luxembourg":               ("Europe", "N"),
    "Malta":                    ("Europe", "N"),
    "Cyprus":                   ("Europe", "N"),

    # Asia
    "China":                    ("Asia", "N"),
    "Japan":                    ("Asia", "N"),
    "Republic of Korea":        ("Asia", "N"),
    "North Korea":              ("Asia", "N"),
    "India":                    ("Asia", "N"),
    "Pakistan":                 ("Asia", "N"),
    "Bangladesh":                ("Asia", "N"),
    "Sri Lanka":                ("Asia", "N"),
    "Nepal":                    ("Asia", "N"),
    "Bhutan":                   ("Asia", "N"),
    "Maldives":                 ("Asia", "N"),
    "Afghanistan":              ("Asia", "N"),
    "Iran":                     ("Asia", "N"),
    "Iraq":                     ("Asia", "N"),
    "Israel":                   ("Asia", "N"),
    "Palestine":                ("Asia", "N"),
    "Jordan":                   ("Asia", "N"),
    "Lebanon":                  ("Asia", "N"),
    "Syrian Arab Republic":     ("Asia", "N"),
    "Saudi Arabia":             ("Asia", "N"),
    "Yemen":                    ("Asia", "N"),
    "United Arab Emirates":     ("Asia", "N"),
    "Oman":                     ("Asia", "N"),
    "Qatar":                    ("Asia", "N"),
    "Bahrain":                  ("Asia", "N"),
    "Kuwait":                   ("Asia", "N"),
    "Turkey":                   ("Asia", "N"),
    "Armenia":                  ("Asia", "N"),
    "Azerbaijan":               ("Asia", "N"),
    "Georgia":                  ("Asia", "N"),
    "Kazakhstan":               ("Asia", "N"),
    "Kyrgyzstan":               ("Asia", "N"),
    "Tajikistan":               ("Asia", "N"),
    "Turkmenistan":             ("Asia", "N"),
    "Uzbekistan":               ("Asia", "N"),
    "Mongolia":                 ("Asia", "N"),
    "Thailand":                 ("Asia", "N"),
    "Vietnam":                  ("Asia", "N"),
    "Cambodia":                 ("Asia", "N"),
    "Laos":                     ("Asia", "N"),
    "Myanmar":                  ("Asia", "N"),
    "Malaysia":                 ("Asia", "N"),
    "Singapore":                ("Asia", "N"),
    "Indonesia":                ("Asia", "S"),  # most population south of equator
    "Philippines":              ("Asia", "N"),
    "Brunei":                   ("Asia", "N"),
    "Timor-Leste":              ("Asia", "S"),
    "Taiwan":                   ("Asia", "N"),
    "Hong Kong":                ("Asia", "N"),
    "Macao":                    ("Asia", "N"),

    # Africa
    "South Africa":             ("Africa", "S"),
    "Egypt":                    ("Africa", "N"),
    "Nigeria":                  ("Africa", "N"),
    "Kenya":                    ("Africa", "S"),  # population centroid south of equator
    "Ethiopia":                 ("Africa", "N"),
    "Algeria":                  ("Africa", "N"),
    "Morocco":                  ("Africa", "N"),
    "Tunisia":                  ("Africa", "N"),
    "Libya":                    ("Africa", "N"),
    "Sudan":                    ("Africa", "N"),
    "South Sudan":              ("Africa", "N"),
    "Senegal":                  ("Africa", "N"),
    "Ghana":                    ("Africa", "N"),
    "Cote d'Ivoire":            ("Africa", "N"),
    "Cameroon":                 ("Africa", "N"),
    "Tanzania":                 ("Africa", "S"),
    "Uganda":                   ("Africa", "N"),
    "Rwanda":                   ("Africa", "S"),
    "Burundi":                  ("Africa", "S"),
    "Mozambique":               ("Africa", "S"),
    "Zambia":                   ("Africa", "S"),
    "Zimbabwe":                 ("Africa", "S"),
    "Angola":                   ("Africa", "S"),
    "Madagascar":               ("Africa", "S"),
    "Mauritius":                ("Africa", "S"),
    "Botswana":                 ("Africa", "S"),
    "Namibia":                  ("Africa", "S"),
    "Eswatini":                 ("Africa", "S"),
    "Lesotho":                  ("Africa", "S"),
    "Mali":                     ("Africa", "N"),
    "Niger":                    ("Africa", "N"),
    "Chad":                     ("Africa", "N"),
    "Central African Republic": ("Africa", "N"),
    "Burkina Faso":             ("Africa", "N"),
    "Guinea":                   ("Africa", "N"),
    "Guinea-Bissau":            ("Africa", "N"),
    "Liberia":                  ("Africa", "N"),
    "Sierra Leone":             ("Africa", "N"),
    "Mauritania":               ("Africa", "N"),
    "Gambia":                   ("Africa", "N"),
    "Benin":                    ("Africa", "N"),
    "Togo":                     ("Africa", "N"),
    "Equatorial Guinea":        ("Africa", "N"),
    "Gabon":                    ("Africa", "S"),
    "Congo":                    ("Africa", "S"),
    "Democratic Republic of the Congo": ("Africa", "S"),
    "Cabo Verde":               ("Africa", "N"),
    "Sao Tome and Principe":    ("Africa", "N"),
    "Comoros":                  ("Africa", "S"),
    "Djibouti":                 ("Africa", "N"),
    "Eritrea":                  ("Africa", "N"),
    "Somalia":                  ("Africa", "N"),
    "Malawi":                   ("Africa", "S"),

    # Oceania
    "Australia":                ("Oceania", "S"),
    "New Zealand":              ("Oceania", "S"),
    "Papua New Guinea":         ("Oceania", "S"),
    "Fiji":                     ("Oceania", "S"),
    "Samoa":                    ("Oceania", "S"),
    "Tonga":                    ("Oceania", "S"),
    "Vanuatu":                  ("Oceania", "S"),
    "Solomon Islands":          ("Oceania", "S"),
    "Kiribati":                 ("Oceania", "N"),
    "Marshall Islands":         ("Oceania", "N"),
    "Micronesia":               ("Oceania", "N"),
    "Palau":                    ("Oceania", "N"),
    "Nauru":                    ("Oceania", "S"),
    "Tuvalu":                   ("Oceania", "S"),
}


# US states (and DC). Used to detect multi-state spread within a single
# country. NCBI's geo_admin1 for US isolates is typically the full state
# name. We accept abbreviations too.
US_STATES = frozenset({
    "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
    "Connecticut", "Delaware", "Florida", "Georgia", "Hawaii", "Idaho",
    "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana",
    "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
    "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada",
    "New Hampshire", "New Jersey", "New Mexico", "New York",
    "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon",
    "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota",
    "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington",
    "West Virginia", "Wisconsin", "Wyoming", "District of Columbia", "DC",
})


def canonical_country(raw: str | None) -> str | None:
    """Normalize a country name to its canonical form. None if blank."""
    if not raw:
        return None
    s = raw.strip()
    if not s:
        return None
    # Try synonym first
    key = s.lower()
    if key in _SYNONYMS:
        return _SYNONYMS[key]
    # Try direct match against canonical names
    if s in _COUNTRY_DATA:
        return s
    # Try case-insensitive direct match
    for canonical in _COUNTRY_DATA.keys():
        if canonical.lower() == key:
            return canonical
    # Unknown country — return as-is, callers can detect via continent_of returning None
    return s


def continent_of(country: str | None) -> str | None:
    """Return continent name for a country, or None if unknown."""
    c = canonical_country(country)
    if c is None:
        return None
    data = _COUNTRY_DATA.get(c)
    return data[0] if data else None


def hemisphere_of(country: str | None) -> str | None:
    """Return 'N' or 'S' for a country, or None if unknown."""
    c = canonical_country(country)
    if c is None:
        return None
    data = _COUNTRY_DATA.get(c)
    return data[1] if data else None


def is_us_state(admin1: str | None) -> bool:
    if not admin1:
        return False
    return admin1.strip() in US_STATES
