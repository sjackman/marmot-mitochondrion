Submit-block ::= {
  contact {
    contact {
      name name {
        last "Jackman",
        first "Shaun",
        initials "S.D.",
      },
      affil std {
        affil "BC Cancer Agency Genome Sciences Centre",
        div "Bioinformatics",
        city "Vancouver",
        sub "BC",
        country "Canada",
        street "100-570 W 7th Ave",
        email "sjackman@gmail.com",
        fax "+1-604-876-3561",
        phone "+1-604-707-5800",
        postal-code "V5Z 4S6"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Jackman",
            first "Shaun",
            initials "S.D.",
          }
        },
        {
          name name {
            last "Birol",
            first "Inanc",
            initials "I.",
          }
        },
        {
          name name {
            last "Jones",
            first "Steven",
            initials "S.J.M.",
          }
        }
      },
      affil std {
        affil "BC Cancer Agency Genome Sciences Centre",
        div "Bioinformatics",
        city "Vancouver",
        sub "BC",
        country "Canada",
        street "100-570 W 7th Ave",
        postal-code "V5Z 4S6"
      }
    }
  },
  subtype new
}

Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Jackman",
              first "Shaun",
              initials "S.D.",
            }
          },
          {
            name name {
              last "Birol",
              first "Inanc",
              initials "I.",
            }
          },
          {
            name name {
              last "Jones",
              first "Steven",
              initials "S.J.M.",
            }
          }
        },
        affil std {
          affil "BC Cancer Agency Genome Sciences Centre",
          div "Bioinformatics",
          city "Vancouver",
          sub "BC",
          country "Canada",
          street "100-570 W 7th Ave",
          postal-code "V5Z 4S6"
        }
      },
      title "Sequencing and assembly of yellow-bellied marmot (Marmota flaviventris)"
    }
  }
}

Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        "PRJNA491472"
      }
    },
    {
      label str "BioSample",
      num 1,
      data strs {
        "SAMN10078668"
      }
    }
  }
}
