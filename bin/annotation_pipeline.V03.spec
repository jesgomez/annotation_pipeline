{
    "jobs": {
        "split_reference*": {
            "time": "00:10:00",
            "priority": "highprio"
        },
        "process_junctions":{
            "time": "03:00:00"
        },
        "augustus_*":{
            "time": "16h",
            "threads": "2"
        },
        "merge*": {
            "time": "00:10:00",
            "priority": "highprio"
        },
        "geneid_*":{
            "time": "3h", 
            "threads": "2"
        },
        "Genemark":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "48h",
            "priority": "lowprio",
            "threads": "8"
        },
        "Genemark-ET":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "48h",
            "priority": "lowprio",
            "threads": "8"
        },
        "Glimmer":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "48h",
            "priority": "lowprio"
        },
        "spaln":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "20h",
            "threads": "8"
        },
        "Clean_transcripts": {
            "time": "00:10:00",
            "priority": "highprio"
        },
        "Pasa":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "72h",
            "priority": "lowprio",
            "threads": "8"
        },
        "Transdecoder":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "12h"
        },
        "EVM_*":{
            "queue": "himem",
            "time": "90h",
            "priority": "lowprio",
            "threads": 16
        },
        "Select_evm":{
            "time": "3h"   
        },
        "Annotation_Update*":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "3h",
            "threads": 2
        },
        "cmsearch":{
            "queue": "himem",
            "time": "48h",
            "priority": "lowprio"
        },
        "ncRNA_*":{
             "time": "00:10:00"
        },
        "Blast*":{
            "time": "3h"
        },
        "lncRNA":{
            "time": "1h"
        }
    } 
}

