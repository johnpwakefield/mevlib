

import pkgutil


def combinelines(limit, prespaces, lines):
    def splitline(line):
        inwords = line.split()
        outwords = []
        while len(inwords):
            nl = []
            for w in inwords:
                if limit - prespaces - 2 - len(" ".join(nl)) - len(w) > 0:
                    nl.append(w)
                else:
                    break
            if not len(nl):
                print(
                    "Cannot break into lines; either lines are too short "
                    "or words are too big."
                )
                raise Exception("Cannot break lines.")
            outwords.append(nl)
            inwords = inwords[len(nl):]
        return [
            (" ".join(nl), i + 1 == len(outwords))
            for i, nl in enumerate(outwords)
        ]
    return "\n".join([
        prespaces * " " + nl + (" &" if not last else "")
        for ol in lines for nl, last in splitline(ol)
    ])


def f03(outfile, matrices, temperatures, verb=False):

    setd  = r"this%d  = {}".format(matrices[0].shape[1])
    setad = r"this%ad = {}".format(matrices[0].shape[0])
    allav = r"allocate(this%axisvals({}))".format(len(temperatures))
    alldv = r"allocate(this%datavals(this%ad, this%d, {}))".format(
        len(temperatures)
    )
    setav = r"this%axisvals = (/ {} /)".format(
        ", ".join(map(str, temperatures))
    )
    setdvs = [
        r"this%datavals(:, :, {}) = reshape((/ {} /), (/ {} /))".format(
            i + 1,
            ", ".join(map(str, matrices[i].flatten('F'))),
            ", ".join(map(str, matrices[i].shape))
        )
        for i in range(len(temperatures))
    ]

    initlines = combinelines(
        79, 8, [setd, setad, allav, alldv, setav] + setdvs
    )

    try:
        tpl = pkgutil.get_data('mevlib', '/templates/mevlookup_template.f03')
        if tpl is None:
            raise FileNotFoundError("template not found")
        program = tpl.decode().replace('!inline_init_here!', initlines)
    except FileNotFoundError:
        print("Could not find template file.")
        exit(1)

    with open(outfile, 'w') as f:
        f.write(program)


def f90(outfile, matrices, temperatures, verb=False):

    setd  = r"d  = {}".format(matrices[0].shape[1])
    setad = r"ad = {}".format(matrices[0].shape[0])
    allav = r"allocate(axisvals({}))".format(len(temperatures))
    alldv = r"allocate(datavals(ad, d, {}))".format(len(temperatures))
    setav = r"axisvals = (/ {} /)".format(", ".join(map(str, temperatures)))
    setdvs = [
        r"datavals(:, :, {}) = reshape((/ {} /), (/ {} /))".format(
            i + 1,
            ", ".join(map(str, matrices[i].flatten('F'))),
            ", ".join(map(str, matrices[i].shape))
        )
        for i in range(len(temperatures))
    ]

    initlines = combinelines(
        79, 8, [setd, setad, allav, alldv, setav] + setdvs
    )

    try:
        tpl = pkgutil.get_data('mevlib', '/templates/mevlookup_template.f03')
        if tpl is None:
            raise FileNotFoundError("template not found")
        program = tpl.decode().replace('!inline_init_here!', initlines)
    except FileNotFoundError:
        print("Could not find template file.")
        exit(1)

    with open(outfile, 'w') as f:
        f.write(program)


