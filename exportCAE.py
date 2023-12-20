def export_to_cae(file_name):
    start_data = """<Indata>
<Id>
<Program>ec701</Program>
<Version>2.1.4</Version>
</Id>"""
    load_start = "\n<Last>\n<Kraft>"
    load_end = "\n</Kraft>\n</Last>"

    Fx = 0;    Fy =339;    Fz = 0;    Mx =150;    My = -150;    Mz = 40
    loads = ''
    for lc_nr in range(5):
        load_temp = f'\n<Idx{lc_nr:03d}>\n<Smxd>{Fx:.1f}</Smxd>\n<Smyd>{Fy:.1f}</Smyd>\n<Smzd>{Fz:.1f}</Smzd>\n<Spxd>{Mx:.1f}</Spxd>\n<Spyd>{My:.1f}</Spyd>\n<Spzd>{Mz:.1f}</Spzd>\n</Idx{lc_nr:03d}>'
        loads+= load_temp
    pile_numb = 8
    pile_start = '\n<Pale>'

    pile_end = '\n</Pale>'

    pile_pos = ''
    for pile_nr in range(10):
        alfa = 180; lengd = 5.1 ; lutnig = 4; typ = 0;xkrd = -0.85;ykrd = 2.85; zkrd = 0;

        pile_pos_temp = f'\n<Idx{pile_nr:03d}>\n<alpha>{alfa:.1f}</alpha>\n<lengd>{lengd:.1f}</lengd>\n<lutning>{lutnig:.1f}</lutning>\n<typ>{typ}</typ>\n<xkrd>{xkrd:.2f}</xkrd>\n<ykrd>{ykrd:.2f}</ykrd>\n<zkrd>{zkrd:.2f}</zkrd>\n</Idx{pile_nr:03d}>'
        pile_pos += pile_pos_temp

    pile_total = f'<antal>{pile_nr+1}</antal>'



    end_data = "\n</Indata>"
    final_file = start_data+load_start+loads+load_end+pile_start+pile_pos+pile_total+pile_end+end_data
    print(final_file)




export_to_cae("xyz.cea")