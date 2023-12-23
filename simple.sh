ualpha:
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o ualpha.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s ualpha

smis_jacobian
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o smis_jacobian_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_jacobian -k 2 -x

smis_all
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o smis_all_k2_stratified.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s smis_all -k 2 -x

cmis
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o cmis.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s cmis

average
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o Average.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s average

proxy
cargo run --features="pbrt progress-bar" --release --example=cli -- -t -2 -n 1 -o proxy.png -m 0.2 "scene/meeting_ply.pbrt" plane-single -n 40960 -s proxy_sample


proxy Elapsed Integrator: 32667 ms
sims Elapsed Integrator: 18328 ms

