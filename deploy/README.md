# Deploying a public BitMako demo

This is the step-by-step for building a small subset of the full 1.36B-compound
corpus and hosting it as a public, clickable demo. Steps 1–2 run on your
workstation (where the real `data/` files already live); steps 3+ run on a VPS.

## 1. Pick a subset size and build it (on your workstation)

20M compounds is a good default: representative, builds in minutes, and the
resulting files (~1 GB index, ~2.5 GB fp-store, ~400 MB prop-store, ~500 MB
skip index) fit comfortably on a small VPS disk. Adjust `--limit` for a
bigger/smaller demo.

```powershell
cargo build --release

# Carve the subset out of the real compounds.lance (no re-ingest needed).
.\target\release\bitmako.exe extract-subset `
  --lance data\compounds.lance --output data\demo\compounds.lance --limit 20000000

# Build the same four artifacts you'd build for the full corpus, just against
# the subset — identical commands, smaller/faster.
.\target\release\bitmako.exe build-index `
  --lance data\demo\compounds.lance --output data\demo\compounds.bitmako
.\target\release\bitmako.exe build-skip `
  --index data\demo\compounds.bitmako --output data\demo\compounds.skip
.\target\release\bitmako.exe build-fp-store `
  --lance data\demo\compounds.lance --output data\demo\compounds.fp
.\target\release\bitmako.exe build-prop-store `
  --lance data\demo\compounds.lance --output data\demo\compounds.prop
```

**Test it locally first** before shipping anything to a server:

```powershell
.\target\release\bitmako.exe serve `
  --index data\demo\compounds.bitmako --skip data\demo\compounds.skip `
  --fp-store data\demo\compounds.fp --prop-store data\demo\compounds.prop `
  --lance data\demo\compounds.lance --port 8080 `
  --demo-notice "Demo instance — a subset of the full 1.36B-compound corpus."
```

Visit `http://localhost:8080/` — confirm search works and the amber demo
banner shows up.

## 2. Pick hosting

Recommended: a small VPS (Hetzner CX22 ~€4/mo, or a DigitalOcean/Linode
$6/mo droplet) running Ubuntu 22.04+. Rather than cross-compiling from
Windows (fragile with Lance's native/protoc dependencies), **build the release
binary on the VPS itself** — a one-time ~15–25 min build, same as on your
workstation.

## 3. Provision the VPS

```bash
ssh root@YOUR_VPS_IP

apt update && apt install -y build-essential pkg-config libssl-dev \
  protobuf-compiler git nginx certbot python3-certbot-nginx

curl https://sh.rustup.rs -sSf | sh -s -- -y
source "$HOME/.cargo/env"

useradd -m -s /usr/sbin/nologin bitmako
mkdir -p /opt/bitmako/{bin,data/demo}
chown -R bitmako:bitmako /opt/bitmako

git clone https://github.com/AndrejRumenovski/bitmako.git ~/bitmako-src
cd ~/bitmako-src
cargo build --release
cp target/release/bitmako /opt/bitmako/bin/
```

## 4. Copy the subset data up

From your workstation:

```powershell
scp -r data\demo\* root@YOUR_VPS_IP:/opt/bitmako/data/demo/
```

Then on the VPS: `chown -R bitmako:bitmako /opt/bitmako/data`

## 5. Install the systemd service

```bash
cp ~/bitmako-src/deploy/bitmako-demo.service /etc/systemd/system/
systemctl daemon-reload
systemctl enable --now bitmako-demo
systemctl status bitmako-demo   # should show "active (running)"
curl http://127.0.0.1:8080/health   # sanity check before exposing publicly
```

## 6. Point DNS and put nginx + TLS in front

Add an `A` record for your chosen subdomain (e.g. `demo.yourdomain.com`)
pointing at the VPS's IP in your domain registrar's DNS panel.

```bash
sed 's/DEMO_DOMAIN/demo.yourdomain.com/' ~/bitmako-src/deploy/nginx-bitmako.conf \
  > /etc/nginx/sites-available/bitmako-demo
ln -s /etc/nginx/sites-available/bitmako-demo /etc/nginx/sites-enabled/
nginx -t && systemctl reload nginx

certbot --nginx -d demo.yourdomain.com   # provisions + auto-renews TLS
```

Visit `https://demo.yourdomain.com/` — you should see the live search UI.

## Note on rate limiting in this topology

BitMako's own per-IP limiter (`src/api.rs`) identifies clients by TCP peer
address. Once nginx is the reverse proxy, every request's peer address is
nginx itself (`127.0.0.1`) — so BitMako's limiter becomes one **shared**
budget across all visitors combined (30 req/min total), not 30-per-visitor.
The nginx config above (`limit_req_zone $binary_remote_addr`) is what actually
enforces per-visitor limits in this deployment; it sees the real client IP.
Both layers still help (nginx = per-visitor, BitMako = a global ceiling), just
know they're not doing the same job.

## Updating the demo later

```bash
cd ~/bitmako-src && git pull && cargo build --release
cp target/release/bitmako /opt/bitmako/bin/
systemctl restart bitmako-demo
```
