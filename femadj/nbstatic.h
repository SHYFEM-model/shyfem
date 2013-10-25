
c how to mark static nodes (not moveable nodes)
c
c	iastatic gives the node type that has to be considered static
c	if you do not want to use static nodes, set it to -1
c	better not use the value 0
c	if set to a positive values it marks all nodes with this
c	node type as static (no exchange, no smooth)
c	boundary nodes are considered static in any case
c
c	nbstatic is the flag used in nbound to indicate static nodes
c	nbstatic cannot be 0 (internal) or 1 (boundary)
c	you can leave the default value below

        integer nbstatic
        parameter (nbstatic = 4)

        integer iastatic
        parameter (iastatic = -1)
        !parameter (iastatic = 1)

