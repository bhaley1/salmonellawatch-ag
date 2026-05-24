"""Simplified continent silhouettes for the cluster map background.

These are deliberately coarse polygon outlines — enough vertices to read
as "world map" at 320×180px, simple enough to stay tiny in the SVG.
Each list is a sequence of (longitude, latitude) tuples forming a closed
ring (the closing vertex is implicit, SVG closepath handles it).

Coordinates approximate the major coastline curves. Inland features,
island chains smaller than ~500km, and minor inlets are dropped. This
is map decoration, not geography.

If anyone needs better-looking continents later, the right move is to
swap this module out for a Natural Earth 1:110m import, but at 320×180px
display size, more detail wouldn't be visible anyway.
"""

from __future__ import annotations


# Each value is a list of (lon, lat) tuples. Some continents have multiple
# polygons (e.g. island chains separated from the main landmass).

NORTH_AMERICA = [
    # Main landmass: Canada + USA + Mexico + Central America
    [
        (-168, 65), (-156, 71), (-128, 71), (-110, 73), (-95, 75), (-80, 74),
        (-65, 60), (-55, 53), (-55, 47), (-64, 45), (-67, 44), (-70, 41),
        (-75, 37), (-80, 32), (-82, 25), (-90, 29), (-93, 30), (-97, 26),
        (-97, 21), (-92, 18), (-87, 16), (-83, 9), (-79, 8), (-82, 12),
        (-94, 17), (-105, 20), (-115, 28), (-118, 33), (-124, 40), (-124, 48),
        (-130, 54), (-140, 60), (-156, 60), (-166, 56), (-168, 65),
    ],
    # Greenland (separate polygon)
    [
        (-55, 60), (-43, 60), (-30, 70), (-22, 73), (-22, 80), (-50, 83),
        (-65, 81), (-55, 78), (-50, 70), (-55, 60),
    ],
    # Cuba + nearby Caribbean simplification
    [
        (-85, 22), (-74, 20), (-75, 23), (-83, 23), (-85, 22),
    ],
]

SOUTH_AMERICA = [
    [
        # Caribbean coast / Colombia
        (-77, 12), (-72, 11), (-62, 11),
        # Venezuela / Guyana coast
        (-52, 5), (-50, 0),
        # Amazon mouth / Brazilian east coast bulge
        (-48, -1), (-44, -3), (-38, -5), (-35, -8), (-37, -12), (-39, -15),
        (-41, -18), (-43, -23), (-48, -25), (-52, -29), (-55, -34),
        (-58, -38), (-66, -45), (-72, -53),
        # Tierra del Fuego
        (-70, -56), (-66, -55),
        # Up the Pacific (Andes) coast
        (-72, -42), (-72, -34), (-71, -28), (-71, -18), (-79, -8),
        (-81, -4), (-78, 2), (-77, 8), (-77, 12),
    ],
]

EUROPE = [
    # Mainland Europe (connected to Asia, but we draw Europe up to the Urals
    # and let Asia start there). Includes Scandinavia and British Isles
    # approximated as the same outline.
    [
        # Start at Gibraltar, go up the Atlantic coast
        (-9, 36), (-9, 43),
        # North coast of Spain to French Atlantic
        (-2, 43), (-1, 46), (-4, 48),
        # Channel and across to British Isles (approximated)
        (-2, 50), (-6, 50), (-10, 55), (-7, 58), (-2, 60),
        # Back over the North Sea
        (3, 56), (5, 58), (8, 63), (12, 67),
        # Northern Scandinavia
        (22, 71), (32, 70), (40, 68),
        # Down to the Urals (border with Asia)
        (52, 70), (66, 68),
        # Eastern border south through Caspian
        (60, 60), (50, 52), (45, 48), (40, 42),
        # Black Sea, Aegean, Mediterranean coast
        (34, 36), (29, 37), (23, 38), (18, 40), (14, 41),
        # Italian boot, French Riviera
        (12, 45), (8, 43), (3, 42),
        # Spanish Mediterranean coast back to Gibraltar
        (-2, 43), (-7, 38), (-9, 36),
    ],
]

AFRICA = [
    [
        (-17, 21), (-17, 14), (-13, 8), (-8, 4), (5, 5), (9, 4), (10, 2),
        (12, -5), (14, -11), (12, -16), (17, -29), (20, -35), (26, -34),
        (32, -28), (35, -22), (40, -16), (41, -10), (51, -1), (43, 11),
        (46, 14), (40, 16), (37, 22), (35, 31), (30, 32), (22, 32), (15, 31),
        (10, 30), (3, 31), (-2, 30), (-9, 28), (-13, 27), (-17, 21),
    ],
]

ASIA = [
    # Mainland Asia: Urals → Pacific → Indian Ocean → Arabian peninsula → back.
    # Walk counterclockwise around the perimeter for one continuous outline.
    [
        # Northern boundary: from the Urals across to the Bering Strait
        (66, 68), (75, 73), (90, 73), (105, 77), (130, 73), (160, 71),
        # East coast: Kamchatka south through Korea
        (170, 62), (155, 58), (140, 54), (132, 47), (122, 40),
        # China east coast
        (123, 33), (122, 25), (115, 22),
        # Indochina peninsula
        (108, 21), (108, 11), (104, 10), (100, 6),
        # Malay peninsula tip and west side
        (100, 1), (103, 1), (98, 8), (98, 12),
        # Bay of Bengal coast → Indian peninsula
        (90, 22), (85, 21), (80, 8), (78, 8), (74, 12),
        # Arabian Sea coast
        (68, 21), (62, 25), (58, 25), (54, 24),
        # Persian Gulf south coast, Arabian Peninsula
        (49, 22), (45, 13), (43, 12), (38, 14),
        # Red Sea east coast, Sinai/Levant
        (36, 18), (35, 27), (35, 30), (36, 31), (35, 35),
        # Turkey, Caucasus, Caspian, central Asia
        (40, 37), (45, 40), (50, 45), (60, 50), (66, 60), (66, 68),
    ],
    # Japan (rough single shape combining Honshu)
    [
        (130, 31), (140, 35), (142, 41), (142, 45), (140, 41), (134, 34),
        (130, 31),
    ],
    # Sri Lanka
    [
        (80, 6), (82, 6), (82, 9), (80, 9), (80, 6),
    ],
    # Borneo
    [
        (109, -3), (118, -3), (119, 4), (109, 7), (109, -3),
    ],
    # Sumatra
    [
        (95, 5), (105, -2), (104, -6), (100, -2), (95, 5),
    ],
    # New Guinea
    [
        (130, -2), (150, -3), (151, -10), (135, -9), (130, -2),
    ],
    # Philippines (rough single shape)
    [
        (118, 5), (122, 6), (125, 12), (122, 18), (120, 16), (118, 10), (118, 5),
    ],
]

AUSTRALIA = [
    [
        (113, -22), (114, -29), (118, -35), (130, -32), (138, -36), (146, -39),
        (149, -37), (153, -28), (153, -25), (146, -20), (142, -11), (135, -12),
        (130, -12), (122, -17), (113, -22),
    ],
    # New Zealand (rough single shape combining both islands)
    [
        (172, -41), (175, -37), (178, -38), (177, -46), (170, -46), (172, -41),
    ],
]


ALL_CONTINENTS = [
    ("north_america", NORTH_AMERICA),
    ("south_america", SOUTH_AMERICA),
    ("europe", EUROPE),
    ("africa", AFRICA),
    ("asia", ASIA),
    ("australia", AUSTRALIA),
]


def all_polygons():
    """Yield (name, polygon) for every continent ring."""
    for name, polys in ALL_CONTINENTS:
        for poly in polys:
            if poly:
                yield (name, poly)
