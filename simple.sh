ualpha:
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 10 -o .\result\ualpha.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 409600 -s ualpha

smis_jacobian
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\smis_jacobian_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_jacobian -k 2 -x

smis_all
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\smis_all_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_all -k 4 -x

cmis
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\cmis.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s cmis

average
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\Average.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s average

proxy
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\proxy.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s proxy_sample

# same sample count
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 3 -o .\result\proxy.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s proxy_sample
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 3 -o .\result\smis_all_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_all -k 2 -x
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o .\result\smis_all_k4_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_all -k 4 -x
$echo