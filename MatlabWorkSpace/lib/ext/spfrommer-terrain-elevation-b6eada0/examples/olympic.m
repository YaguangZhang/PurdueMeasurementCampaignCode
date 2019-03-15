addpath('../')

region = fetchregion([47.328158, 48.436983], [-124.983574, -122.852227], ...
                     'display', true);
elevData = region.readelevation([47.328158, 48.436983], ...
                                [-124.983574, -122.852227], ...
                                'SampleFactor', 100, ...
                                'display', true);
dispelev(elevData, 'mode', 'cartesain');