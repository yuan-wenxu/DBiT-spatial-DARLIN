#!/bin/bash

chip_preset_names() {
    printf '%s\n' "50-50 50-20 100-20"
}

chip_preset_names_csv() {
    printf '%s\n' "50-50, 50-20, 100-20"
}

chip_preset_is_supported() {
    case "$1" in
        50-50|50-20|100-20) return 0 ;;
        *) return 1 ;;
    esac
}

apply_chip_preset() {
    local selected_chip=$1
    case "$selected_chip" in
        50-50)
            chip=50-50
            x_spots_number=50
            y_spots_number=50
            length_spot=50
            interval=50
            whitelist_path="$REPO_DIR/docs/barcodes/barcodes.tsv"

            ;;
        50-20)
            chip=50-20
            x_spots_number=50
            y_spots_number=50
            length_spot=20
            interval=20
            whitelist_path="$REPO_DIR/docs/barcodes/barcodes.tsv"

            ;;
        100-20)
            chip=100-20
            x_spots_number=100
            y_spots_number=100
            length_spot=20
            interval=20
            whitelist_path="$REPO_DIR/docs/barcodes/barcodes100.tsv"

            ;;
        *)
            return 1
            ;;
    esac

}
