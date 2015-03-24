{
    "jobs": {
        "split_reference": {
            "time": "00:10:00"
        },
        "process_junctions":{
            "time": "03:00:00"
        },
        "augustus_*":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "16h"
        },
        "merge*": {
            "time": "00:10:00"
        },
        "geneid_*":{
            "time": "3h"
        },
        "Genemark":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "12h"
        },
        "Genemark-ET":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "12h"
        },
        "Glimmer":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "12h"
        },
        "spaln":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "6h"
        },
        "Clean_transcripts": {
            "time": "00:10:00"
        },
        "Pasa":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "48h",
            "priority": "lowprio"
        },
        "Transdecoder":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "3h"
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
        "Annotation_Update":{
            "queue": "development",
            "account": "AA",
            "extra": ["-C","AA"],
            "time": "10h"
        }
    } 
}

